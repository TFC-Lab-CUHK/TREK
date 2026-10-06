#!/usr/bin/env python3
"""
GTF Processor: Parse GTF files and build transcript structures for polyA identification
"""

import logging
from typing import Dict, List, Tuple, Set
from collections import defaultdict
from dataclasses import dataclass
from interlap import InterLap

logger = logging.getLogger(__name__)


@dataclass
class Exon:
    """Represents an exon with 1-based genomic coordinates"""
    start: int  # 1-based inclusive
    end: int    # 1-based inclusive
    
    def __post_init__(self):
        if self.start > self.end:
            raise ValueError(f"Invalid exon: start ({self.start}) > end ({self.end})")


@dataclass
class Transcript:
    """Represents a transcript with exon structure"""
    transcript_id: str
    gene_name: str
    chromosome: str
    strand: str
    exons: List[Exon]
    ncbi_gene_id: str = ''
    transcript_biotype: str = ''
    reference_gene_id: str = ''
    
    def __post_init__(self):
        # Sort exons by genomic position
        self.exons.sort(key=lambda e: e.start)
        
    @property
    def is_multi_exon(self) -> bool:
        """Check if transcript has multiple exons"""
        return len(self.exons) > 1

    @property
    def gene_id(self) -> str:
        """NCBI GeneID when available, otherwise the annotation's gene ID."""
        return self.ncbi_gene_id or self.reference_gene_id
    
    @property
    def start(self) -> int:
        """Transcript start position (1-based)"""
        return self.exons[0].start if self.exons else 0
    
    @property
    def end(self) -> int:
        """Transcript end position (1-based)"""
        return self.exons[-1].end if self.exons else 0
    
    def get_splice_junctions(self) -> Tuple[int, ...]:
        """
        Get splice junction coordinates (exon boundaries excluding transcript ends)
        Returns tuple of (end1, start2, end2, start3, ...) for multi-exon transcripts
        Returns empty tuple for single-exon transcripts
        """
        if len(self.exons) <= 1:
            return tuple()
        
        junctions = []
        for i in range(len(self.exons) - 1):
            junctions.append(self.exons[i].end)      # Donor site
            junctions.append(self.exons[i + 1].start)  # Acceptor site
        
        return tuple(junctions)


class GTFProcessor:
    """Process GTF file to extract transcript structures"""
    
    def __init__(self):
        """
        Initialize GTF processor
        """
        self.stats = defaultdict(int)
    
    def parse_gtf(self, gtf_file: str) -> Dict[str, Transcript]:
        """
        Parse GTF file and build transcript structures
        
        Args:
            gtf_file: Path to GTF file
            
        Returns:
            Dictionary of transcript_id -> Transcript
        """
        logger.info(f"Parsing GTF file: {gtf_file}")
        
        # Temporary storage
        transcript_metadata = {}
        transcript_exons = defaultdict(list)
        
        with open(gtf_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                # Skip comments and empty lines
                if line.startswith('#') or not line.strip():
                    continue
                
                self.stats['total_lines'] += 1
                
                try:
                    parsed = self._parse_gtf_line(line.strip())
                    if not parsed:
                        continue
                    
                    feature_type = parsed['feature']
                    
                    # Only process transcript and exon features
                    if feature_type not in {'transcript', 'exon', 'mRNA'}:
                        continue
                    
                    transcript_id = parsed['attributes'].get('transcript_id')
                    if not transcript_id:
                        self.stats['missing_transcript_id'] += 1
                        continue
                    
                    # Handle transcript/mRNA features
                    if feature_type in {'transcript', 'mRNA'}:
                        attrs = parsed['attributes']
                        previous = transcript_metadata.get(transcript_id, {})
                        transcript_metadata[transcript_id] = {
                            'chromosome': parsed['seqname'],
                            'strand': parsed['strand'],
                            # gene_id is used when gene_name is absent
                            'gene_name': attrs.get('gene_name') or previous.get('gene_name') or attrs.get('gene_id', ''),
                            'ncbi_gene_id': attrs.get('ncbi_gene_id') or previous.get('ncbi_gene_id', ''),
                            'reference_gene_id': attrs.get('gene_id') or previous.get('reference_gene_id', ''),
                            'transcript_biotype': attrs.get('transcript_biotype') or previous.get('transcript_biotype', ''),
                        }
                        self.stats['transcripts_found'] += 1
                    
                    # Handle exon features
                    elif feature_type == 'exon':
                        exon = Exon(parsed['start'], parsed['end'])
                        transcript_exons[transcript_id].append(exon)
                        
                        # Extract metadata from exon if transcript not seen
                        if transcript_id not in transcript_metadata:
                            attrs = parsed['attributes']
                            transcript_metadata[transcript_id] = {
                                'chromosome': parsed['seqname'],
                                'strand': parsed['strand'],
                                'gene_name': attrs.get('gene_name') or attrs.get('gene_id', ''),
                                'ncbi_gene_id': attrs.get('ncbi_gene_id', ''),
                                'reference_gene_id': attrs.get('gene_id', ''),
                                'transcript_biotype': attrs.get('transcript_biotype', ''),
                            }
                        else:
                            # Some annotations supply gene metadata only on exons.
                            attrs = parsed['attributes']
                            metadata = transcript_metadata[transcript_id]
                            for target, source in (('reference_gene_id', 'gene_id'),
                                                   ('ncbi_gene_id', 'ncbi_gene_id'),
                                                   ('gene_name', 'gene_name'),
                                                   ('transcript_biotype', 'transcript_biotype')):
                                if not metadata.get(target) and attrs.get(source):
                                    metadata[target] = attrs[source]
                
                except Exception as e:
                    logger.warning(f"Line {line_num}: Error parsing - {e}")
                    self.stats['parsing_errors'] += 1
        
        # Build transcript objects
        transcripts = {}
        for transcript_id, metadata in transcript_metadata.items():
            exons = transcript_exons.get(transcript_id, [])
            
            if not exons:
                self.stats['transcripts_no_exons'] += 1
                continue
            
            try:
                transcript = Transcript(
                    transcript_id=transcript_id,
                    gene_name=metadata['gene_name'],
                    chromosome=metadata['chromosome'],
                    strand=metadata['strand'],
                    exons=exons,
                    ncbi_gene_id=metadata.get('ncbi_gene_id', ''),
                    transcript_biotype=metadata.get('transcript_biotype', ''),
                    reference_gene_id=metadata.get('reference_gene_id', ''),
                )
                transcripts[transcript_id] = transcript
                self.stats['valid_transcripts'] += 1
                
                if transcript.is_multi_exon:
                    self.stats['multi_exon_transcripts'] += 1
                else:
                    self.stats['single_exon_transcripts'] += 1
                    
            except Exception as e:
                logger.error(f"Error creating transcript {transcript_id}: {e}")
                self.stats['transcript_errors'] += 1
        
        self._log_stats()
        return transcripts
    
    def _parse_gtf_line(self, line: str) -> Dict:
        """Parse a single GTF line"""
        fields = line.split('\t')
        if len(fields) != 9:
            return None
        
        # Parse attributes column
        attributes = {}
        for attr in fields[8].split(';'):
            attr = attr.strip()
            if not attr:
                continue
            
            if ' "' in attr and attr.endswith('"'):
                key, value = attr.split(' "', 1)
                key = key.strip()
                value = value.rstrip('"')
                # db_xref appears multiple times; capture the GeneID entry specifically
                if key == 'db_xref':
                    if value.startswith('GeneID:'):
                        attributes['ncbi_gene_id'] = value[len('GeneID:'):]
                else:
                    attributes[key] = value
        
        return {
            'seqname': fields[0],
            'feature': fields[2],
            'start': int(fields[3]),
            'end': int(fields[4]),
            'strand': fields[6],
            'attributes': attributes
        }
    
    def organize_transcripts(self, transcripts: Dict[str, Transcript]) -> Tuple[
        Dict[str, Dict[Tuple[int, ...], Tuple[str, int, int, str]]],  # multi_exon_dict
        Dict[Tuple[str, str], InterLap]                           # single_exon_dict
    ]:
        """
        Organize transcripts for efficient read assignment
        
        Returns:
            - multi_exon_dict: Nested dict mapping chr -> sj_pattern -> (transcript_id, start, end, strand)
            - single_exon_dict: InterLap objects per (chr, strand) containing (start, end, transcript_id)
        """
        logger.info("Organizing transcripts for read assignment")
        
        multi_exon_dict = defaultdict(dict)
        single_exon_dict = defaultdict(InterLap)
        
        for transcript_id, transcript in transcripts.items():
            chrom_strand_key = (transcript.chromosome, transcript.strand)
            
            if transcript.is_multi_exon:
                # Map splice junction pattern to transcript (single value), keyed by chrom only
                sj_pattern = transcript.get_splice_junctions()
                if sj_pattern in multi_exon_dict[transcript.chromosome]:
                    logger.warning(f"Duplicate splice junction pattern for {transcript_id} and {multi_exon_dict[transcript.chromosome][sj_pattern][0]}")
                multi_exon_dict[transcript.chromosome][sj_pattern] = (
                    transcript_id,
                    transcript.start,
                    transcript.end,
                    transcript.strand
                )
            else:
                # Single exon transcript - add to InterLap
                exon = transcript.exons[0]
                single_exon_dict[chrom_strand_key].add((exon.start, exon.end, transcript_id))
        
        logger.info(f"Organized {sum(len(v) for v in multi_exon_dict.values())} unique splice patterns across {len(multi_exon_dict)} chromosomes")
        logger.info(f"Organized {len(single_exon_dict)} InterLap objects for single-exon transcripts")
        
        return dict(multi_exon_dict), dict(single_exon_dict)
    
    def _log_stats(self):
        """Log processing statistics"""
        logger.info("GTF Processing Statistics:")
        logger.info(f"  Total lines processed: {self.stats['total_lines']}")
        logger.info(f"  Valid transcripts: {self.stats['valid_transcripts']}")
        logger.info(f"    Multi-exon: {self.stats['multi_exon_transcripts']}")
        logger.info(f"    Single-exon: {self.stats['single_exon_transcripts']}")
        
        if self.stats['parsing_errors'] > 0:
            logger.warning(f"  Parsing errors: {self.stats['parsing_errors']}")
        if self.stats['transcripts_no_exons'] > 0:
            logger.warning(f"  Transcripts without exons: {self.stats['transcripts_no_exons']}")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
