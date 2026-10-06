#!/usr/bin/env python3
"""
Alignment Processor: Handle minimap2 alignment and read-to-transcript assignment
"""

import os
import subprocess
import tempfile
import logging
import pysam
from interlap import InterLap
from typing import Dict, List, Tuple, Optional
from collections import defaultdict
from tqdm import tqdm

logger = logging.getLogger(__name__)

PLATFORM_PROFILES = {
    'drna': {'preset': 'splice', 'splice_strand': 'f', 'kmer': 14},
    'pcr-cdna': {'preset': 'splice', 'splice_strand': 'b', 'kmer': 15},
    'dcdna': {'preset': 'splice', 'splice_strand': 'b', 'kmer': 15},
    # Read orientation is not assumed for cDNA or PacBio libraries;
    # --splice-strand f declares reads to be in transcript orientation.
    'pacbio': {'preset': 'splice:hq', 'splice_strand': 'b', 'kmer': 15},
}
DEFAULT_ALIGNMENT_PROFILE = {'preset': 'splice', 'splice_strand': 'b', 'kmer': 15}
UNSPECIFIED_PLATFORMS = {'', '-', '.', 'na', 'n/a', 'unknown', 'not reported', 'unspecified'}


def normalize_platform(value):
    """Normalise a library-type label; empty or NA-like values mean unspecified."""
    if value is None:
        return 'unspecified'
    if not isinstance(value, str):
        raise ValueError('Platform must be a library-type name or empty')
    label = value.strip().lower()
    if label in UNSPECIFIED_PLATFORMS:
        return 'unspecified'
    if label not in PLATFORM_PROFILES:
        raise ValueError(f'Unknown platform {value!r}; use drna, pcr-cdna, dcdna, pacbio, or leave it empty')
    return label


def alignment_profiles(fastq_files, platform=None, splice_strand=None):
    """Resolve per-file types; missing types use splice -u b -k15."""
    supplied = [platform] if isinstance(platform, str) or platform is None else list(platform)
    supplied = [normalize_platform(value) for value in (supplied or [None])]
    if len(supplied) == 1:
        supplied *= len(fastq_files)
    elif len(supplied) != len(fastq_files):
        raise ValueError('Provide one library type for all FASTQs, or one per FASTQ in the same order')
    if isinstance(splice_strand, str):
        splice_strand = splice_strand.lower()
    if splice_strand not in (None, 'f', 'b'):
        raise ValueError('splice_strand must be f or b')
    resolved = []
    for path, label in zip(fastq_files, supplied):
        settings = dict(DEFAULT_ALIGNMENT_PROFILE if label == 'unspecified' else PLATFORM_PROFILES[label])
        if splice_strand is not None:
            settings['splice_strand'] = splice_strand
        resolved.append(dict(fastq=os.path.basename(str(path)), platform=label, **settings, max_intron=500000))
    return resolved


# pysam CIGAR operation codes
_OP_M = 0   # Match/mismatch
_OP_D = 2   # Deletion
_OP_N = 3   # Intron (splice junction)
_OP_EQ = 7  # Sequence match
_OP_X = 8   # Sequence mismatch


def _extract_junctions(cigartuples, ref_start: int) -> Optional[Tuple[int, ...]]:
    """
    Extract splice junctions from pysam cigartuples.
    
    Args:
        cigartuples: List of (op, length) from pysam
        ref_start: 0-based reference start
        
    Returns:
        Tuple of junction coordinates (1-based) or None for single-block reads
    """
    junctions = []
    pos = ref_start  # 0-based
    for op, length in cigartuples:
        if op == _OP_M or op == _OP_EQ or op == _OP_X:
            pos += length
        elif op == _OP_N:
            junctions.append(pos)          # Donor: 0-based end = 1-based last position
            pos += length
            junctions.append(pos + 1)      # Acceptor: 1-based first position
        elif op == _OP_D:
            pos += length
        # I(1), S(4), H(5), P(6) don't advance reference position
    return tuple(junctions) if junctions else None


_4GB = 4 * 1024 ** 3  # 4 GiB in bytes


def _get_genome_size(fasta_path: str) -> int:
    """
    Return the total genome size in bytes by summing sequence lengths in a FASTA file.
    Only sequence characters (non-header lines) are counted.
    """
    total = 0
    with open(fasta_path, 'r') as fh:
        for line in fh:
            if not line.startswith('>'):
                total += len(line.rstrip('\n\r'))
    return total


class AlignmentProcessor:
    """Process alignments and assign reads to transcripts"""
    
    def __init__(self,
                 multi_exon_dict: Dict,
                 single_exon_dict: Dict,
                 min_mapq: int = 1,
                 min_overlap_ratio: float = 0.5,
                 terminal_exon_dict: Optional[Dict] = None):
        """
        Initialize alignment processor
        
        Args:
            multi_exon_dict: Multi-exon transcript dictionary from GTFProcessor
            single_exon_dict: Single exon dictionary from GTFProcessor
            min_mapq: Minimum mapping quality
            min_overlap_ratio: Minimum overlap ratio (intersection/transcript_length) for single-exon assignment
            terminal_exon_dict: Terminal exon intervals for every annotated
                multi-exon model, keyed by (chromosome, annotation strand).
        """
        if not 0 < min_overlap_ratio <= 1:
            raise ValueError('min_overlap_ratio must be greater than 0 and at most 1')
        self.multi_exon_dict = multi_exon_dict
        self.single_exon_dict = single_exon_dict
        if terminal_exon_dict is None:
            # Build terminal-exon intervals from the chain index when none are supplied.
            terminal_exon_dict = defaultdict(InterLap)
            for chrom, models in multi_exon_dict.items():
                for chain, (tid, start, end, strand) in models.items():
                    terminal_start, terminal_end = (chain[-1], end) if strand == '+' else (start, chain[0])
                    terminal_exon_dict[chrom, strand].add((terminal_start, terminal_end, tid))
        self.terminal_exon_dict = terminal_exon_dict
        self.min_mapq = min_mapq
        self.min_overlap_ratio = min_overlap_ratio
        self.stats = defaultdict(int)
    
    def process_alignment(self,
                         genome_fasta: str,
                         fastq_files: List[str],
                         bed_file: str,
                         threads: int = 8,
                         platform=None,
                         splice_strand: Optional[str] = None) -> Dict[str, List[int]]:
        """
        Run minimap2 and process alignments on-the-fly.
        Select minimap2 settings from explicit sequencing/library types.
        Types can be supplied once for all files, or separately for each file.
        
        Args:
            genome_fasta: Path to genome FASTA
            fastq_files: List of FASTQ files
            bed_file: Path to reference BED file (converted from GTF)
            threads: Number of threads
            platform: Per-file library types; missing values use the generic RNA profile
            splice_strand: Override the library defaults with minimap2 -u f/b
            
        Returns:
            Dictionary mapping transcript_id to list of read_end_positions
        """
        profiles = alignment_profiles(fastq_files, platform, splice_strand)
        # Genomes above 4 GiB are indexed in one chunk by raising -I above the genome size.
        genome_size = _get_genome_size(genome_fasta)
        split_args: List[str] = []
        if genome_size > _4GB:
            index_size_gb = int(genome_size / 1024 ** 3) + 4  # genome size + 4 GiB headroom
            split_args = ['-I', f'{index_size_gb}g']
            logger.info(
                f"Genome size {genome_size / 1024**3:.2f} GiB > 4 GiB: "
                f"adding -I{index_size_gb}g to minimap2 to avoid index splitting"
            )
        else:
            logger.info(f"Genome size {genome_size / 1024**3:.2f} GiB: default -I is sufficient")

        transcript_reads: Dict[str, list] = defaultdict(list)
        groups = {}
        for path, profile in zip(fastq_files, profiles):
            key = (profile['platform'], profile['preset'], profile['splice_strand'], profile['kmer'])
            groups.setdefault(key, []).append(str(path))
        for (label, preset, direction, kmer), files in groups.items():
            logger.info(f"Running minimap2 for {label}: {len(files)} file(s), "
                        f"preset={preset}, -u {direction}, -k {kmer}")
            self._run_minimap2(
                genome_fasta=genome_fasta,
                fastq_files=files,
                bed_file=bed_file,
                threads=threads,
                preset=preset,
                extra_args=['-u', direction, '-k', str(kmer), '-G', '500000'] + split_args,
                desc=f"Processing {label} reads",
                transcript_reads=transcript_reads,
                single_exon_oriented=(direction == 'f'),
            )
        
        self._log_stats()
        logger.info(f"Assigned reads to {len(transcript_reads)} transcripts")
        
        # Sort read ends within each transcript.
        sorted_transcript_reads = {
            transcript_id: sorted(positions)
            for transcript_id, positions in transcript_reads.items()
        }
        
        return sorted_transcript_reads

    def _run_minimap2(self,
                      genome_fasta: str,
                      fastq_files: List[str],
                      bed_file: str,
                      threads: int,
                      preset: str,
                      extra_args: List[str],
                      desc: str,
                      transcript_reads: Dict,
                      single_exon_oriented: bool = False) -> None:
        """
        Run a single minimap2 invocation and process SAM output on-the-fly,
        accumulating results into the shared transcript_reads dict.
        
        Args:
            genome_fasta: Path to genome FASTA
            fastq_files: FASTQ files to align
            bed_file: BED junction guide file (may be empty string / None)
            threads: Number of threads for minimap2
            preset: minimap2 preset string (e.g. 'splice' or 'splice:hq')
            extra_args: Additional minimap2 flags inserted before reference/query
            desc: Label for the tqdm progress bar
            transcript_reads: Shared dict that results are written into
            single_exon_oriented: Whether input reads are known to follow RNA
                orientation (the effective -u f profile).
        """
        base_cmd = ['minimap2', '-ax', preset, '-t', str(threads)] + extra_args + ['--secondary=no']
        
        if bed_file:
            base_cmd.extend(['--junc-bed', bed_file])
            logger.info(f"Using junction guide from: {bed_file}")
        
        base_cmd.append(genome_fasta)
        
        # Assignments from a file are merged after its aligner exits successfully.
        for file_idx, fastq_file in enumerate(fastq_files, 1):
            cmd = base_cmd + [fastq_file]
            file_desc = f"{desc} [{file_idx}/{len(fastq_files)}] {fastq_file}"
            logger.info(f"Command: {' '.join(cmd)}")
            
            file_assignments = defaultdict(list)
            # Aligner stderr is kept in a temporary file; its tail is reported on failure.
            with tempfile.TemporaryFile() as stderr_log:
                process = None
                try:
                    process = subprocess.Popen(
                        cmd, stdout=subprocess.PIPE, stderr=stderr_log, text=False,
                    )
                    with pysam.AlignmentFile(process.stdout, 'r') as sam:
                        for read in tqdm(sam, desc=file_desc):
                            # The counters below are mutually exclusive and sum to alignment_records.
                            self.stats['alignment_records'] += 1

                            # Unmapped, secondary and supplementary records
                            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                                self.stats['filtered_alignments'] += 1
                                continue

                            # Mapped primary records below the mapping-quality threshold
                            if read.mapping_quality < self.min_mapq:
                                self.stats['low_mapq'] += 1
                                continue

                            # Primary segments of split alignments (SA tag)
                            if read.has_tag('SA'):
                                self.stats['sa_tagged_primary'] += 1
                                continue

                            # Mapped primary records that enter transcript assignment
                            self.stats['total_reads_processed'] += 1
                            
                            chrom = read.reference_name
                            ref_start = read.reference_start   # 0-based
                            ref_end = read.reference_end         # 0-based
                            strand = '-' if read.is_reverse else '+'
                            
                            junctions = _extract_junctions(read.cigartuples, ref_start)
                            assigned = self._assign_read(
                                chrom, strand, ref_start, ref_end, junctions,
                                single_exon_oriented=single_exon_oriented,
                            )
                            
                            if not assigned:
                                self.stats['unassigned_reads'] += 1
                                continue
                            
                            assigned_transcript, transcript_strand = assigned
                            
                            # Determine 3' end position (1-based) based on transcript strand
                            if transcript_strand == '+':
                                read_end = ref_end
                            else:
                                read_end = ref_start + 1
                            
                            file_assignments[assigned_transcript].append(read_end)
                            self.stats['assigned_reads'] += 1
                    process.wait()
                    if process.returncode != 0:
                        raise RuntimeError(f'minimap2 exited with status {process.returncode}')
                except BaseException as error:
                    if process is not None:
                        if process.poll() is None:
                            process.kill()
                        process.wait()
                    self.stats['failed_files'] += 1
                    if not isinstance(error, Exception):
                        raise
                    stderr_log.seek(0, os.SEEK_END)
                    stderr_log.seek(max(0, stderr_log.tell() - 8192))
                    detail = stderr_log.read().decode(errors='replace').strip()
                    raise RuntimeError(
                        f'Mapping failed for {fastq_file}: {error}' + (f'\n{detail}' if detail else '')
                    ) from error
                finally:
                    if process is not None and process.stdout is not None:
                        process.stdout.close()
            for tid, positions in file_assignments.items():
                transcript_reads[tid].extend(positions)

    def _assign_read(self, chrom: str, strand: str,
                     ref_start: int, ref_end: int,
                     junctions: Optional[Tuple[int, ...]],
                     single_exon_oriented: bool = False) -> Optional[Tuple[str, str]]:
        """
        Assign a single read to a transcript
        
        Args:
            chrom: Chromosome name
            strand: Read strand ('+' or '-')
            ref_start: 0-based reference start
            ref_end: 0-based reference end
            junctions: Splice junction tuple from _extract_junctions, or None
            single_exon_oriented: Restrict unspliced reads to the alignment strand
                only when the input follows the transcript's RNA orientation.
        
        Returns:
            Tuple of (transcript_id, transcript_strand), or None if unassigned
        """
        if junctions:
            # Multi-exon read - match based on splice pattern (strand-agnostic)
            chrom_multi_exon_dict = self.multi_exon_dict.get(chrom)
            if chrom_multi_exon_dict:
                match = chrom_multi_exon_dict.get(junctions)
                if match:
                    transcript_id, tx_start, tx_end, tx_strand = match
                    self.stats['multi_exon_assigned'] += 1
                    return transcript_id, tx_strand
                else:
                    self.stats['multi_exon_no_match'] += 1
            else:
                self.stats['multi_exon_no_match'] += 1
        else:
            # Unspliced reads: candidate strands follow the library orientation;
            # a read eligible on both strands is left unassigned.
            candidate_strands = (strand,) if single_exon_oriented else ('+', '-')
            read_start, read_end = ref_start + 1, ref_end  # 1-based inclusive
            eligible = []
            for tx_strand in candidate_strands:
                interlap = self.single_exon_dict.get((chrom, tx_strand))
                if interlap is None:
                    continue
                best_transcript = None
                best_ratio = 0.0
                best_intersection, best_length = 0, 0
                for tx_start, tx_end, transcript_id in interlap.find((read_start, read_end)):
                    intersect_length = max(0, min(read_end, tx_end) - max(read_start, tx_start) + 1)
                    transcript_length = tx_end - tx_start + 1
                    overlap_ratio = intersect_length / transcript_length
                    # Highest fraction wins within a strand; ties keep the first candidate.
                    if overlap_ratio > best_ratio:
                        best_ratio = overlap_ratio
                        best_transcript = transcript_id
                        best_intersection, best_length = intersect_length, transcript_length
                if best_transcript is not None and best_ratio >= self.min_overlap_ratio:
                    eligible.append((best_transcript, tx_strand, best_intersection, best_length))

            if len(eligible) == 1:
                tid, assigned_strand, single_overlap, single_length = eligible[0]
                # Exclude the read if any multi-exon terminal exon overlaps it by a
                # fraction of its own length >= the single-exon candidate's fraction.
                for tx_strand in candidate_strands:
                    terminal_exons = self.terminal_exon_dict.get((chrom, tx_strand))
                    if terminal_exons is None:
                        continue
                    for exon_start, exon_end, _ in terminal_exons.find((read_start, read_end)):
                        exon_overlap = max(0, min(read_end, exon_end) - max(read_start, exon_start) + 1)
                        exon_length = exon_end - exon_start + 1
                        if exon_overlap * single_length >= single_overlap * exon_length:
                            self.stats['single_exon_terminal_overlap'] += 1
                            return None
                self.stats['single_exon_assigned'] += 1
                return tid, assigned_strand
            if len(eligible) > 1:
                self.stats['single_exon_ambiguous_strand'] += 1
            else:
                self.stats['single_exon_no_overlap'] += 1
        
        return None
    
    def _log_stats(self):
        """Log processing statistics"""
        logger.info("Read Assignment Statistics:")
        logger.info(f"  Alignment records read: {self.stats['alignment_records']}")
        logger.info(f"  Excluded before assignment:")
        logger.info(f"    Unmapped, secondary or supplementary: {self.stats['filtered_alignments']}")
        logger.info(f"    Low MAPQ (mapped primary): {self.stats['low_mapq']}")
        logger.info(f"    SA-tagged primary alignments: {self.stats['sa_tagged_primary']}")
        logger.info(f"  Reads entering transcript assignment: {self.stats['total_reads_processed']}")
        logger.info(f"    Assigned reads: {self.stats['assigned_reads']}")
        logger.info(f"      Multi-exon: {self.stats['multi_exon_assigned']}")
        logger.info(f"      Single-exon: {self.stats['single_exon_assigned']}")
        logger.info(f"    Unassigned reads: {self.stats['unassigned_reads']}")
        logger.info(f"      Single-exon reads with eligible candidates on both strands: "
                    f"{self.stats['single_exon_ambiguous_strand']}")
        logger.info(f"      Single-exon reads with equal or better terminal-exon overlap fractions: "
                    f"{self.stats['single_exon_terminal_overlap']}")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
