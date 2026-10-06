#!/usr/bin/env python3
"""
TREK: per-sample pipeline for identifying polyadenylation sites from long-read data

This pipeline processes long-read sequencing data to identify alternative 
transcription end sites (TES) / polyadenylation sites:
1. Process GTF annotation to extract transcript structures
2. Align reads using minimap2 with splice-aware mapping
3. Assign reads to transcripts based on splice junction patterns
4. Identify alternative TES using Gaussian Mixture Models
"""

import argparse
import csv
import json
import logging
import sys
import os
import tempfile
import pickle
from pathlib import Path
from collections import defaultdict
from contextlib import contextmanager
from interlap import InterLap

if __package__:
    from .gtf_processor import GTFProcessor
    from .gff_processor import GFFProcessor
    from .alignment_processor import AlignmentProcessor, alignment_profiles, normalize_platform
    from .apa_finder import TESFinder
    from .internal_priming_filter import InternalPrimingFilter
else:
    from gtf_processor import GTFProcessor
    from gff_processor import GFFProcessor
    from alignment_processor import AlignmentProcessor, alignment_profiles, normalize_platform
    from apa_finder import TESFinder
    from internal_priming_filter import InternalPrimingFilter

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)
ASSIGNMENT_FILTER_SCHEMA = 'trek.read_assignment_filters.v1'


def published_file_mode():
    """Regular-file permissions under the current umask."""
    current = os.umask(0)
    os.umask(current)
    return 0o666 & ~current


@contextmanager
def atomic_output(path, mode='w'):
    """Publish a file only after its contents have been written successfully."""
    path = Path(path)
    temporary = tempfile.NamedTemporaryFile(
        mode=mode, encoding=None if 'b' in mode else 'utf-8',
        dir=path.parent, prefix=f'.{path.name}.', suffix='.tmp', delete=False,
    )
    try:
        with temporary as handle:
            yield handle
        os.chmod(temporary.name, published_file_mode())
        os.replace(temporary.name, path)
    finally:
        Path(temporary.name).unlink(missing_ok=True)


def read_input_list(filename):
    """Read FASTQ paths and optional types; relative paths use the list directory."""
    manifest = Path(filename).expanduser().resolve()
    fastq_files, platforms, seen = [], [], set()
    first_row = True
    with manifest.open(encoding='utf-8-sig', newline='') as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip('\r\n')
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            location = f'{manifest}:{line_number}'
            try:
                # Double quotes group a path containing spaces; backslashes are literal.
                if '\t' in line:
                    fields = next(csv.reader([line], delimiter='\t', skipinitialspace=True, strict=True))
                else:
                    fields = [field for field in next(csv.reader(
                        [line.strip()], delimiter=' ', skipinitialspace=True, strict=True)) if field != '']
            except (ValueError, csv.Error) as error:
                raise ValueError(f'{location}: malformed input-list row: {error}') from error
            if len(fields) not in (1, 2):
                raise ValueError(f'{location}: expected FASTQ path and optional platform (one or two columns)')
            fields = [field.strip() for field in fields]
            if (first_row and len(fields) == 2 and fields[0].lower() in ('fastq', 'fastq_path')
                    and fields[1].lower() in ('platform', 'library_type')):
                first_row = False
                continue
            first_row = False
            if not fields[0]:
                raise ValueError(f'{location}: the FASTQ path is empty')
            try:
                platform = normalize_platform(fields[1] if len(fields) == 2 else None)
            except ValueError as error:
                raise ValueError(f'{location}: {error}') from error
            fastq = Path(fields[0]).expanduser()
            if not fastq.is_absolute():
                fastq = manifest.parent/fastq
            fastq = fastq.resolve()
            if not fastq.is_file():
                raise FileNotFoundError(f'{location}: FASTQ file not found: {fastq}')
            if fastq in seen:
                raise ValueError(f'{location}: duplicate FASTQ path: {fastq}')
            seen.add(fastq)
            fastq_files.append(str(fastq))
            platforms.append(platform)
    if not fastq_files:
        raise ValueError(f'{manifest}: no FASTQ entries found')
    return fastq_files, platforms


class TrekPipeline:
    """Per-sample pipeline: annotation, alignment, assignment, PAS calling, priming filter."""
    
    def __init__(self, gtf_file, genome_fasta, fastq_files, output_dir, prefix='polyA',
                 threads=8, min_mapq=1, min_overlap_ratio=0.5,
                 min_reads=10, min_cluster_size=10, max_clusters=5, min_distance=50,
                 min_relative_dominance=0.3, min_sharpness=0.5, n_jobs=-1,
                 priming_window=20, priming_a_threshold=0.5,
                 random_seed=42, filter_priming=True, platform=None, splice_strand=None):
        """Initialize pipeline with parameters"""
        self.gtf_file = gtf_file
        self.genome_fasta = genome_fasta
        self.fastq_files = fastq_files
        self.platform = platform
        self.splice_strand = splice_strand
        self.alignment_profiles = alignment_profiles(fastq_files, platform, splice_strand)
        # Inputs are identified by their resolved paths.
        self.input_identity = {
            'annotation': str(Path(gtf_file).expanduser().resolve()),
            'genome_fasta': str(Path(genome_fasta).expanduser().resolve()),
            'fastq_files': [str(Path(path).expanduser().resolve()) for path in fastq_files],
        }
        self.output_dir = Path(output_dir)
        self.prefix = prefix
        self.threads = threads
        self.min_mapq = min_mapq
        self.min_overlap_ratio = min_overlap_ratio
        if not 0 < min_overlap_ratio <= 1:
            raise ValueError('min_overlap_ratio must be greater than 0 and at most 1')
        self.assignment_policy = {
            'single_exon_strand_rule': 'same_strand_if_u_f_else_unique_eligible_strand_v1',
            'single_exon_terminal_exon_filter': 'terminal_fraction_ge_best_single_fraction_v1',
            'min_overlap_ratio': min_overlap_ratio,
            'min_mapq': min_mapq,
        }
        self.min_reads = min_reads
        self.min_cluster_size = min_cluster_size
        self.max_clusters = max_clusters
        self.min_distance = min_distance
        self.min_relative_dominance = min_relative_dominance
        self.min_sharpness = min_sharpness
        self.n_jobs = n_jobs
        self.priming_window = priming_window
        self.priming_a_threshold = priming_a_threshold
        self.filter_priming = filter_priming
        self.random_seed = random_seed
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info("TREK pipeline started")
        logger.info(f"Output directory: {self.output_dir}")
    
    def run(self):
        """Run the complete pipeline"""
        try:
            logger.info("STEP 1: Processing GTF annotation")
            transcripts, multi_exon_dict, single_exon_dict = self._process_gtf()

            assignment_file = self.output_dir / f"{self.prefix}.read_assignments.pkl"
            transcript_reads = self._load_assignments_if_valid(assignment_file)

            if transcript_reads is not None:
                logger.info(f"Found intact assignment file, skipping alignment (Steps 2-3): {assignment_file}")
            else:
                logger.info("STEP 2: Converting GTF to BED for alignment guide")
                bed_file = self._gtf_to_bed(transcripts)

                logger.info("STEP 3: Running alignment and assigning reads to transcripts")
                transcript_reads = self._process_alignment(bed_file, multi_exon_dict, single_exon_dict, transcripts)
            
            logger.info("STEP 4: Identifying alternative polyA sites")
            apa_results = self._find_apa_sites(transcript_reads)

            removed_apa = {}
            if self.filter_priming:
                logger.info("STEP 5: Filtering internal priming artifacts")
                apa_results, removed_apa = self._filter_internal_priming(
                    apa_results, transcripts
                )
            else:
                logger.info("STEP 5: Internal priming filtering disabled; skipping")
            
            logger.info("STEP 6: Writing results")
            self._write_results(apa_results, transcripts)
            self._write_removed_apa_results(removed_apa, transcripts)
            
            logger.info("Pipeline completed successfully!")
            
        except Exception as e:
            logger.error(f"Pipeline failed: {e}", exc_info=True)
            raise
    
    def _load_assignments_if_valid(self, assignment_file):
        """Reuse assignments only with matching filtering and assignment rules."""
        if not assignment_file.exists():
            return None
        filter_file = assignment_file.with_suffix('.filters.json')
        cache_error = (
            f"Cannot reuse assignment cache without confirmed SA-tag filtering: {assignment_file}. "
            "Cached read-end positions do not contain SA tags. Use a new output directory "
            "or prefix to regenerate assignments with SA-tagged reads excluded."
        )
        try:
            filters = json.loads(filter_file.read_text())
        except (OSError, ValueError) as error:
            raise ValueError(cache_error) from error
        if (not isinstance(filters, dict)
                or filters.get('schema') != ASSIGNMENT_FILTER_SCHEMA
                or filters.get('exclude_sa_tagged_reads') is not True):
            raise ValueError(cache_error)
        if filters.get('alignment_profiles') != self.alignment_profiles:
            raise ValueError(
                f"Cannot reuse assignment cache with different or unrecorded platform settings: {assignment_file}. "
                "Use a new output directory or prefix for the requested library types and splice-strand settings."
            )
        if filters.get('assignment_policy') != self.assignment_policy:
            raise ValueError(
                f"Cannot reuse assignment cache with different or unrecorded assignment rules: {assignment_file}. "
                "Use a new output directory or prefix to apply the current single-exon strand, overlap and terminal-exon rules."
            )
        if filters.get('mapping_completed') is not True:
            raise ValueError(
                f'Cannot reuse assignment cache without confirmed mapping completion: {assignment_file}. '
                'Use a new output directory or prefix to regenerate complete assignments.'
            )
        if filters.get('input_identity') != self.input_identity:
            raise ValueError(
                f'Cannot reuse assignment cache with different or unrecorded input/reference paths: {assignment_file}. '
                'Use a new output directory or prefix for these inputs and references.'
            )
        try:
            with open(assignment_file, 'rb') as f:
                data = pickle.load(f)
            logger.info(f"Loaded {len(data)} transcript assignments from cache")
            return data
        except Exception as e:
            logger.warning(f"Assignment file exists but could not be loaded ({e}), re-running alignment")
            return None

    def _process_gtf(self):
        """Process GTF or GFF3 file"""
        logger.info(f"Processing annotation: {self.gtf_file}")
        
        ext = Path(self.gtf_file).suffix.lower()
        if ext in ('.gff', '.gff3'):
            processor = GFFProcessor()
            transcripts = processor.parse_gff(self.gtf_file)
        else:
            processor = GTFProcessor()
            transcripts = processor.parse_gtf(self.gtf_file)
        
        logger.info(f"Parsed {len(transcripts)} transcripts")
        if not transcripts:
            raise ValueError(
                f'Annotation yielded no transcript models: {self.gtf_file}. '
                'Check the file format (GTF needs transcript/exon rows with transcript_id; '
                'GFF3 needs exon rows with Parent) and that it is not empty.'
            )
        
        multi_exon_dict, single_exon_dict = processor.organize_transcripts(transcripts)
        return transcripts, multi_exon_dict, single_exon_dict
    
    def _gtf_to_bed(self, transcripts=None):
        """Build BED12 from the same parsed GTF/GFF models used for assignment."""
        if transcripts is None:
            transcripts, _, _ = self._process_gtf()
        if not transcripts:
            return None
        bed_file = self.output_dir / f"{self.prefix}.junctions.bed"
        with atomic_output(bed_file) as handle:
            for tx in sorted(transcripts.values(), key=lambda t: (t.chromosome, t.start, t.transcript_id)):
                start = tx.start - 1
                sizes = ','.join(str(e.end - e.start + 1) for e in tx.exons) + ','
                offsets = ','.join(str(e.start - 1 - start) for e in tx.exons) + ','
                fields = (tx.chromosome, start, tx.end, tx.transcript_id, 0, tx.strand,
                          start, tx.end, 0, len(tx.exons), sizes, offsets)
                handle.write('\t'.join(map(str, fields)) + '\n')
        logger.info(f"Created BED file: {bed_file}")
        return str(bed_file)
    
    def _process_alignment(self, bed_file, multi_exon_dict, single_exon_dict, transcripts):
        """Run alignment and assign reads on-the-fly"""
        # Terminal exons of every annotated multi-exon model.
        terminal_exon_dict = defaultdict(InterLap)
        for tid, transcript in transcripts.items():
            if transcript.is_multi_exon:
                exon = transcript.exons[-1] if transcript.strand == '+' else transcript.exons[0]
                terminal_exon_dict[transcript.chromosome, transcript.strand].add((exon.start, exon.end, tid))
        processor = AlignmentProcessor(
            multi_exon_dict=multi_exon_dict,
            single_exon_dict=single_exon_dict,
            min_mapq=self.min_mapq,
            min_overlap_ratio=self.min_overlap_ratio,
            terminal_exon_dict=dict(terminal_exon_dict),
        )
        
        transcript_reads = processor.process_alignment(
            genome_fasta=self.genome_fasta,
            fastq_files=self.fastq_files,
            bed_file=bed_file,
            threads=self.threads,
            platform=self.platform,
            splice_strand=self.splice_strand,
        )
        
        # Save read assignments as pickle
        assignment_file = self.output_dir / f"{self.prefix}.read_assignments.pkl"
        with atomic_output(assignment_file, 'wb') as f:
            pickle.dump(transcript_reads, f)
        metadata = {
            'schema': ASSIGNMENT_FILTER_SCHEMA,
            'exclude_sa_tagged_reads': True,
            'mapping_completed': True,
            'input_identity': self.input_identity,
            'alignment_profiles': self.alignment_profiles,
            'assignment_policy': self.assignment_policy,
            'sa_tagged_primary_alignments_removed': processor.stats.get('sa_tagged_primary', 0),
            'single_exon_ambiguous_strand_reads_removed': processor.stats.get('single_exon_ambiguous_strand', 0),
            'single_exon_terminal_overlap_reads_removed': processor.stats.get('single_exon_terminal_overlap', 0),
            'filter_stage': 'After MAPQ and alignment-flag filtering, before transcript assignment',
        }
        with atomic_output(assignment_file.with_suffix('.filters.json')) as handle:
            json.dump(metadata, handle, indent=2)
            handle.write('\n')
        
        logger.info(f"Saved read assignments: {assignment_file}")
        return transcript_reads
    
    def _find_apa_sites(self, transcript_reads):
        """Find alternative polyA sites"""
        finder = TESFinder(
            min_reads=self.min_reads,
            min_cluster_size=self.min_cluster_size,
            max_k=self.max_clusters,
            min_distance=self.min_distance,
            min_relative_dominance=self.min_relative_dominance,
            min_sharpness=self.min_sharpness,
            n_jobs=self.n_jobs,
            random_seed=self.random_seed
        )
        
        return finder.find_apa_sites(transcript_reads)
    
    def _filter_internal_priming(self, apa_results, transcripts):
        """Filter APA sites to remove internal priming artifacts"""
        filter_obj = InternalPrimingFilter(
            genome_fasta=self.genome_fasta,
            window_size=self.priming_window,
            a_content_threshold=self.priming_a_threshold
        )
        
        return filter_obj.filter_apa_results(apa_results, transcripts)

    def _write_removed_apa_results(self, removed_apa, transcripts):
        """Write APA sites removed by the internal priming filter."""
        output_file = self.output_dir / f"{self.prefix}.internal_priming_removed.txt"

        with open(output_file, 'w') as f:
            f.write(
                "transcript_id\tgene_id\tgene_name\tchromosome\tstrand\t"
                "ID\tsite_position\tsite_count\tsite_abundance\t"
                "transcript_biotype\ta_content\n"
            )

            for transcript_id, apa in removed_apa.items():
                transcript = transcripts.get(transcript_id)
                if not transcript:
                    continue

                for position, count, abundance, a_content in zip(
                    apa.site, apa.count, apa.abundance, apa.a_content
                ):
                    locus_id = (
                        f"{transcript.chromosome}:{position}:{transcript.strand}"
                    )
                    f.write(
                        f"{transcript_id}\t{transcript.gene_id}\t"
                        f"{transcript.gene_name}\t{transcript.chromosome}\t"
                        f"{transcript.strand}\t{locus_id}\t{position}\t{count}\t"
                        f"{abundance:.4f}\t{transcript.transcript_biotype}\t"
                        f"{a_content:.4f}\n"
                    )

        logger.info(f"Saved removed internal priming sites: {output_file}")
    
    def _write_results(self, apa_results, transcripts):
        """Write results to output files"""
        # Main results file
        results_file = self.output_dir / f"{self.prefix}.apa_sites.txt"
        
        with open(results_file, 'w') as f:
            f.write("transcript_id\tgene_id\tgene_name\tchromosome\tstrand\t"
                   "ID\tsite_position\tsite_count\tsite_abundance\ttranscript_biotype\n")
            
            for transcript_id, apa in apa_results.items():
                transcript = transcripts.get(transcript_id)
                if not transcript:
                    continue
                
                # Write one line per APA site
                for position, count, abundance in zip(apa.site, apa.count, apa.abundance):
                    # Create locus ID: chrom:position:strand
                    locus_id = f"{transcript.chromosome}:{position}:{transcript.strand}"
                    
                    f.write(f"{transcript_id}\t{transcript.gene_id}\t{transcript.gene_name}\t"
                           f"{transcript.chromosome}\t{transcript.strand}\t{locus_id}\t"
                           f"{position}\t{count}\t{abundance:.4f}\t{transcript.transcript_biotype}\n")
        
        logger.info(f"Saved results: {results_file}")
        
        # Summary file
        summary_file = self.output_dir / f"{self.prefix}.summary.txt"
        
        with open(summary_file, 'w') as f:
            total = len(apa_results)
            with_apa = sum(1 for apa in apa_results.values() if len(apa.site) > 1)
            
            f.write("TREK Summary\n")
            f.write("=" * 50 + "\n")
            f.write(f"Total transcripts analyzed: {total}\n")
            f.write(f"Transcripts with alternative TES: {with_apa}\n")
            percentage = 100 * with_apa / total if total else 0.0
            f.write(f"Percentage with APA: {percentage:.2f}%\n")
            
            site_counts = {}
            for apa in apa_results.values():
                n = len(apa.site)
                site_counts[n] = site_counts.get(n, 0) + 1
            
            f.write("\nDistribution of TES per transcript:\n")
            for n in sorted(site_counts.keys()):
                f.write(f"  {n} sites: {site_counts[n]} transcripts\n")
        
        logger.info(f"Saved summary: {summary_file}")


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='TREK: identify polyadenylation sites from long-read RNA-seq, per transcript'
    )
    
    # Required arguments
    parser.add_argument('-g', '--gtf', required=True, 
                       help='Reference annotation file (GTF, GFF, or GFF3)')
    parser.add_argument('-f', '--fasta', required=True, 
                       help='Genome FASTA file')
    parser.add_argument('-i', '--input', required=True,
                       help='Input list: FASTQ path and optional platform (drna, pcr-cdna, dcdna, pacbio); '
                            'missing platform uses splice -u b -k15')
    
    # Output arguments
    parser.add_argument('-o', '--output', default='results', 
                       help='Output directory (default: results)')
    parser.add_argument('-p', '--prefix', default='polyA', 
                       help='Output file prefix (default: polyA)')
    
    # Alignment arguments
    parser.add_argument('-t', '--threads', type=int, default=8, 
                       help='Number of threads (default: 8)')
    parser.add_argument('--min-mapq', type=int, default=1, 
                       help='Minimum MAPQ (default: 1)')
    parser.add_argument('--min-overlap', type=float, default=0.5,
                       help='Minimum fraction of an annotated single-exon transcript covered (default: 0.5)')
    parser.add_argument('--splice-strand', choices=('f', 'b'), type=str.lower,
                       help='Override minimap2 -u for all inputs: f for known transcript orientation, b for both; '
                            'default: f for drna, b for all other or unspecified types')
    
    # TES detection arguments
    parser.add_argument('--min-reads', type=int, default=10, 
                       help='Minimum reads per transcript (default: 10)')
    parser.add_argument('--min-cluster-size', type=int, default=10, 
                       help='Minimum cluster size (default: 10)')
    parser.add_argument('--max-clusters', type=int, default=5, 
                       help='Maximum clusters (default: 5)')
    parser.add_argument('--min-distance', type=int, default=50, 
                       help='Minimum distance between sites in bp (default: 50)')
    parser.add_argument('--min-dominance', type=float, default=0.3, 
                       help='Minimum relative dominance (default: 0.3)')
    parser.add_argument('--min-sharpness', type=float, default=0.5, 
                       help='Minimum peak sharpness (default: 0.5)')
    
    # Performance arguments
    parser.add_argument('-j', '--jobs', type=int, default=-1, 
                       help='Parallel jobs (default: -1, all CPUs)')
    
    # Internal priming filter arguments
    parser.add_argument('--priming-window', type=int, default=20,
                       help='Downstream transcript-direction window size for A-content check (default: 20)')
    parser.add_argument('--priming-a-threshold', type=float, default=0.5,
                       help='Maximum A-content threshold (default: 0.5)')
    parser.add_argument('--no-filter-priming', action='store_true',
                       help='Disable internal priming filtering (enabled by default)')
    
    # Reproducibility arguments
    parser.add_argument('--serial', action='store_true',
                       help='Run in serial mode for full reproducibility (overrides -j)')
    
    return parser.parse_args()


def main():
    """Main entry point"""
    args = parse_arguments()
    
    try:
        # Validate inputs
        if not Path(args.gtf).exists():
            raise FileNotFoundError(f"Annotation file not found: {args.gtf}")
        if not Path(args.fasta).exists():
            raise FileNotFoundError(f"FASTA file not found: {args.fasta}")
        fastq_files, platforms = read_input_list(args.input)
        
        # Run pipeline
        pipeline = TrekPipeline(
            gtf_file=args.gtf,
            genome_fasta=args.fasta,
            fastq_files=fastq_files,
            platform=platforms,
            splice_strand=args.splice_strand,
            output_dir=args.output,
            prefix=args.prefix,
            threads=args.threads,
            min_mapq=args.min_mapq,
            min_overlap_ratio=args.min_overlap,
            min_reads=args.min_reads,
            min_cluster_size=args.min_cluster_size,
            max_clusters=args.max_clusters,
            min_distance=args.min_distance,
            min_relative_dominance=args.min_dominance,
            min_sharpness=args.min_sharpness,
            n_jobs=1 if args.serial else args.jobs,
            priming_window=args.priming_window,
            priming_a_threshold=args.priming_a_threshold,
            filter_priming=not args.no_filter_priming,
            random_seed=42
        )
        pipeline.run()
        
        return 0
        
    except Exception as e:
        logger.error(f"TREK failed: {e}")
        return 1


if __name__ == '__main__':
    sys.exit(main())
