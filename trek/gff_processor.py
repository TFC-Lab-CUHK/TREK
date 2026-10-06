#!/usr/bin/env python3
"""
GFF/GFF3 Processor: Parse GFF3 files and build transcript structures for polyA identification

GFF3 format uses key=value attributes (vs GTF's key "value"), with hierarchical
relationships expressed through ID/Parent attributes. This processor handles the
gene -> mRNA/transcript -> exon hierarchy to produce the same output structures
as GTFProcessor.
"""

import logging
from urllib.parse import unquote
from typing import Dict, List, Tuple, Set, Optional
from collections import defaultdict
from interlap import InterLap

if __package__:
    from .gtf_processor import Transcript, Exon
else:
    from gtf_processor import Transcript, Exon

logger = logging.getLogger(__name__)

# Feature types recognized as gene-level entries
_GENE_FEATURES = {'gene', 'pseudogene'}

# Features that may carry an ID but are never transcripts. Any other
# ID-bearing, non-gene feature that is the Parent of an exon is a transcript
# (mRNA, lnc_RNA, snoRNA, tRNA, rRNA, miRNA, V_gene_segment, ...).
_NON_TRANSCRIPT_FEATURES = {
    'exon', 'cds', 'intron',
    'five_prime_utr', 'three_prime_utr', 'start_codon', 'stop_codon',
    'region', 'match', 'cdna_match', 'protein_match', 'nucleotide_match',
    'sequence_feature', 'repeat_region', 'origin_of_replication',
}


class GFFProcessor:
    """Process GFF/GFF3 file to extract transcript structures"""

    def __init__(self):
        self.stats = defaultdict(int)

    def parse_gff(self, gff_file: str) -> Dict[str, Transcript]:
        """
        Parse a GFF3 file and build transcript structures.

        Args:
            gff_file: Path to GFF3 file

        Returns:
            Dictionary of transcript_id -> Transcript
        """
        logger.info(f"Parsing GFF3 file: {gff_file}")

        # Storage keyed by GFF3 ID attribute
        gene_info: Dict[str, dict] = {}            # gene ID -> {name, ncbi_gene_id, biotype}
        transcript_metadata: Dict[str, dict] = {}   # transcript ID -> metadata dict
        transcript_exons: Dict[str, List[Exon]] = defaultdict(list)
        exons_by_parent: Dict[str, List[Exon]] = defaultdict(list)

        with open(gff_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#') or not line.strip():
                    continue

                self.stats['total_lines'] += 1

                try:
                    parsed = self._parse_gff_line(line.strip())
                    if not parsed:
                        continue

                    feature = parsed['feature']
                    attrs = parsed['attributes']
                    gff_id = attrs.get('ID', '')
                    parents = attrs.get('Parent', '').split(',') if attrs.get('Parent') else []

                    # ---- gene-level features ----
                    if feature.lower() in {g.lower() for g in _GENE_FEATURES} or feature in _GENE_FEATURES:
                        if gff_id:
                            gene_info[gff_id] = {
                                'name': attrs.get('Name', gff_id),
                                'reference_gene_id': attrs.get('gene_id') or gff_id,
                                'ncbi_gene_id': self._extract_ncbi_gene_id(attrs),
                                'biotype': attrs.get('biotype', attrs.get('gene_biotype', '')),
                                'chromosome': parsed['seqname'],
                                'strand': parsed['strand'],
                            }
                        self.stats['genes_found'] += 1
                        continue

                    # ---- exon-level features ----
                    if feature.lower() == 'exon':
                        exon = Exon(parsed['start'], parsed['end'])

                        if not parents:
                            self.stats['exons_no_parent'] += 1
                            continue

                        for parent in parents:
                            exons_by_parent[parent.strip()].append(exon)
                        continue

                    # ---- candidate transcript-level features ----
                    # Candidates become transcripts once exons are attached to them.
                    if feature.lower() in _NON_TRANSCRIPT_FEATURES:
                        continue
                    tid = self._clean_id(gff_id) if gff_id else None
                    if not tid:
                        self.stats['missing_transcript_id'] += 1
                        continue

                    transcript_metadata[tid] = {
                        'gff_id': gff_id,  # original un-cleaned ID for parent matching
                        'feature': feature,
                        'parents': parents,
                        'attributes': attrs,
                        'chromosome': parsed['seqname'],
                        'strand': parsed['strand'],
                    }

                except Exception as e:
                    logger.warning(f"Line {line_num}: Error parsing - {e}")
                    self.stats['parsing_errors'] += 1

        # Resolve parent relationships after the whole file has been read.
        for parent, exons in exons_by_parent.items():
            tid = self._find_transcript_for_parent(parent, transcript_metadata)
            if tid:
                transcript_exons[tid].extend(exons)
            else:
                self.stats['exons_unresolved_parent'] += len(exons)

        # Keep exon-bearing candidates only.
        transcript_metadata = {
            tid: metadata for tid, metadata in transcript_metadata.items()
            if tid in transcript_exons
        }
        self.stats['transcripts_found'] = len(transcript_metadata)

        # Gene attributes are inherited through nested parents
        # (e.g. miRNA -> primary_transcript -> gene).
        parent_lookup = {metadata['gff_id']: tid for tid, metadata in transcript_metadata.items()}
        # Resolve all gene lineages first; parent lists are consumed below.
        lineage = {
            tid: self._resolve_gene_lineage(metadata['parents'], gene_info, transcript_metadata, parent_lookup)
            for tid, metadata in transcript_metadata.items()
        }
        for tid, metadata in transcript_metadata.items():
            attrs = metadata.pop('attributes')
            parents = metadata.pop('parents')
            metadata.pop('feature')
            parent_gene = lineage[tid]
            gene_parent_id = next((p.strip() for p in parents if p.strip() in gene_info), None)
            metadata['gene_name'] = (
                (parent_gene.get('name', '') if parent_gene else '')
                or attrs.get('Name', attrs.get('gene', metadata['gff_id']))
            )
            metadata['ncbi_gene_id'] = (
                self._extract_ncbi_gene_id(attrs)
                or (parent_gene.get('ncbi_gene_id', '') if parent_gene else '')
            )
            metadata['transcript_biotype'] = (
                attrs.get('biotype', attrs.get('transcript_biotype', ''))
                or (parent_gene.get('biotype', '') if parent_gene else '')
            )
            metadata['reference_gene_id'] = (
                attrs.get('gene_id')
                or (parent_gene.get('reference_gene_id') if parent_gene else None)
                or gene_parent_id
                or (parents[0].strip() if len(parents) == 1 else '')
            )

        # Build Transcript objects
        transcripts = self._build_transcripts(transcript_metadata, transcript_exons)

        self._log_stats()
        return transcripts

    # ------------------------------------------------------------------
    # Parsing helpers
    # ------------------------------------------------------------------

    def _parse_gff_line(self, line: str) -> Optional[Dict]:
        """Parse a single GFF3 line."""
        fields = line.split('\t')
        if len(fields) != 9:
            return None

        attributes = self._parse_gff_attributes(fields[8])

        return {
            'seqname': fields[0],
            'feature': fields[2],
            'start': int(fields[3]),
            'end': int(fields[4]),
            'strand': fields[6],
            'attributes': attributes,
        }

    @staticmethod
    def _parse_gff_attributes(attr_string: str) -> Dict[str, str]:
        """
        Parse GFF3 attribute column (key=value pairs separated by ';').
        Values are URL-decoded per the GFF3 specification.
        """
        attributes: Dict[str, str] = {}
        for item in attr_string.split(';'):
            item = item.strip()
            if not item or '=' not in item:
                continue
            key, _, value = item.partition('=')
            attributes[key.strip()] = unquote(value.strip())
        return attributes

    @staticmethod
    def _clean_id(raw_id: str) -> str:
        """
        Remove common GFF3 ID prefixes such as 'transcript:', 'gene:', 'rna-', etc.
        """
        for prefix in ('transcript:', 'rna-', 'rna:', 'mrna:', 'mRNA:'):
            if raw_id.startswith(prefix):
                return raw_id[len(prefix):]
        return raw_id

    @staticmethod
    def _extract_ncbi_gene_id(attrs: Dict[str, str]) -> str:
        """
        Extract NCBI Gene ID from GFF3 attributes.
        Checks Dbxref / db_xref for 'GeneID:' entries as well as explicit
        'ncbi_gene_id' key. Generic gene IDs are retained separately.
        """
        # Direct attribute
        if 'ncbi_gene_id' in attrs:
            return attrs['ncbi_gene_id']

        # Dbxref (GFF3 standard) or db_xref
        for key in ('Dbxref', 'db_xref', 'dbxref'):
            if key in attrs:
                for ref in attrs[key].split(','):
                    ref = ref.strip()
                    if ref.startswith('GeneID:'):
                        return ref[len('GeneID:'):]
        return ''

    def _resolve_gene_parent(
        self, parents: List[str], gene_info: Dict[str, dict]
    ) -> Optional[dict]:
        """Return gene_info dict for the first matching parent, or None."""
        for p in parents:
            p = p.strip()
            if p in gene_info:
                return gene_info[p]
        return None

    def _resolve_gene_lineage(
        self,
        parents: List[str],
        gene_info: Dict[str, dict],
        transcript_metadata: Dict[str, dict],
        parent_lookup: Dict[str, str],
        max_depth: int = 5,
    ) -> Optional[dict]:
        """Return gene_info for the nearest gene ancestor, walking through
        intermediate transcript features (e.g. miRNA -> primary_transcript -> gene)."""
        frontier = [p.strip() for p in parents if p.strip()]
        seen = set()
        for _ in range(max_depth):
            if not frontier:
                return None
            gene = self._resolve_gene_parent(frontier, gene_info)
            if gene is not None:
                return gene
            next_frontier = []
            for parent in frontier:
                if parent in seen:
                    continue
                seen.add(parent)
                tid = parent_lookup.get(parent) or (
                    self._clean_id(parent) if self._clean_id(parent) in transcript_metadata else None
                )
                if tid is not None:
                    next_frontier.extend(
                        p.strip() for p in transcript_metadata[tid].get('parents', []) if p.strip()
                    )
            frontier = next_frontier
        return None

    def _find_transcript_for_parent(
        self, parent: str, transcript_metadata: Dict[str, dict]
    ) -> Optional[str]:
        """
        Given a Parent value from an exon, find the matching transcript tid.
        Matches against both the cleaned tid and the original gff_id.
        """
        cleaned = self._clean_id(parent)
        if cleaned in transcript_metadata:
            return cleaned
        # Search by original gff_id
        for tid, meta in transcript_metadata.items():
            if meta.get('gff_id') == parent:
                return tid
        return None

    def _build_transcripts(
        self,
        transcript_metadata: Dict[str, dict],
        transcript_exons: Dict[str, List[Exon]],
    ) -> Dict[str, Transcript]:
        """Build Transcript objects from parsed metadata and exons."""
        transcripts: Dict[str, Transcript] = {}

        for tid, metadata in transcript_metadata.items():
            exons = transcript_exons.get(tid, [])
            if not exons:
                self.stats['transcripts_no_exons'] += 1
                continue

            try:
                transcript = Transcript(
                    transcript_id=tid,
                    gene_name=metadata['gene_name'],
                    chromosome=metadata['chromosome'],
                    strand=metadata['strand'],
                    exons=exons,
                    ncbi_gene_id=metadata.get('ncbi_gene_id', ''),
                    transcript_biotype=metadata.get('transcript_biotype', ''),
                    reference_gene_id=metadata.get('reference_gene_id', ''),
                )
                transcripts[tid] = transcript
                self.stats['valid_transcripts'] += 1

                if transcript.is_multi_exon:
                    self.stats['multi_exon_transcripts'] += 1
                else:
                    self.stats['single_exon_transcripts'] += 1

            except Exception as e:
                logger.error(f"Error creating transcript {tid}: {e}")
                self.stats['transcript_errors'] += 1

        return transcripts

    def organize_transcripts(self, transcripts: Dict[str, Transcript]) -> Tuple[
        Dict[str, Dict[Tuple[int, ...], Tuple[str, int, int, str]]],
        Dict[Tuple[str, str], InterLap],
    ]:
        """
        Organize transcripts for efficient read assignment.
        Identical interface to GTFProcessor.organize_transcripts().

        Returns:
            - multi_exon_dict: chr -> sj_pattern -> (transcript_id, start, end, strand)
            - single_exon_dict: (chr, strand) -> InterLap with (start, end, transcript_id)
        """
        logger.info("Organizing transcripts for read assignment")

        multi_exon_dict = defaultdict(dict)
        single_exon_dict = defaultdict(InterLap)

        for transcript_id, transcript in transcripts.items():
            chrom_strand_key = (transcript.chromosome, transcript.strand)

            if transcript.is_multi_exon:
                sj_pattern = transcript.get_splice_junctions()
                if sj_pattern in multi_exon_dict[transcript.chromosome]:
                    logger.warning(
                        f"Duplicate splice junction pattern for {transcript_id} "
                        f"and {multi_exon_dict[transcript.chromosome][sj_pattern][0]}"
                    )
                multi_exon_dict[transcript.chromosome][sj_pattern] = (
                    transcript_id,
                    transcript.start,
                    transcript.end,
                    transcript.strand,
                )
            else:
                exon = transcript.exons[0]
                single_exon_dict[chrom_strand_key].add(
                    (exon.start, exon.end, transcript_id)
                )

        logger.info(
            f"Organized {sum(len(v) for v in multi_exon_dict.values())} "
            f"unique splice patterns across {len(multi_exon_dict)} chromosomes"
        )
        logger.info(
            f"Organized {len(single_exon_dict)} InterLap objects "
            f"for single-exon transcripts"
        )

        return dict(multi_exon_dict), dict(single_exon_dict)

    def _log_stats(self):
        """Log processing statistics."""
        logger.info("GFF3 Processing Statistics:")
        logger.info(f"  Total lines processed: {self.stats['total_lines']}")
        logger.info(f"  Genes found: {self.stats['genes_found']}")
        logger.info(f"  Valid transcripts: {self.stats['valid_transcripts']}")
        logger.info(f"    Multi-exon: {self.stats['multi_exon_transcripts']}")
        logger.info(f"    Single-exon: {self.stats['single_exon_transcripts']}")

        if self.stats['parsing_errors'] > 0:
            logger.warning(f"  Parsing errors: {self.stats['parsing_errors']}")
        if self.stats['transcripts_no_exons'] > 0:
            logger.warning(
                f"  Transcripts without exons: {self.stats['transcripts_no_exons']}"
            )



if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
    )
