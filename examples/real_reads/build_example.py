#!/usr/bin/env python3
"""Rebuild the small real-read example from an indexed BAM and full references."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess

import pysam

from trek.alignment_processor import _extract_junctions
from trek.gtf_processor import GTFProcessor

TRANSCRIPTS = ('NM_181697.3', 'NM_001289746.2', 'NM_002046.7')
CHROMOSOMES = {'NM_181697.3': 'NC_000001.11',
               'NM_001289746.2': 'NC_000012.12', 'NM_002046.7': 'NC_000012.12'}
ROOT = Path(__file__).resolve().parents[2]


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def reference_records(gtf):
    parser = GTFProcessor()
    selected = []
    chromosomes = tuple((chrom + '\t').encode() for chrom in sorted(set(CHROMOSOMES.values())))
    tokens = [f'transcript_id "{tid}"'.encode() for tid in TRANSCRIPTS]
    digest = hashlib.sha256()
    with gtf.open('rb') as handle:
        for raw in handle:
            digest.update(raw)
            if not raw.startswith(chromosomes) or not any(token in raw for token in tokens):
                continue
            line = raw.decode().rstrip('\n\r')
            record = parser._parse_gtf_line(line)
            tid = record['attributes'].get('transcript_id') if record else None
            if tid in CHROMOSOMES and record['seqname'] == CHROMOSOMES[tid]:
                selected.append((tid, line.split('\t'), record))
    if {tid for tid, _, _ in selected} != set(TRANSCRIPTS):
        raise ValueError('Reference does not contain the expected three transcript models')
    return selected, digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bam', type=Path, required=True)
    parser.add_argument('--fasta', type=Path, required=True)
    parser.add_argument('--gtf', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--reads-per-transcript', type=int, default=100)
    parser.add_argument('--seed', type=int, default=20261004)
    parser.add_argument('--flank', type=int, default=2000)
    parser.add_argument('--min-mapq', type=int, default=20)
    args = parser.parse_args()
    if args.reads_per_transcript < 10 or args.flank < 100:
        raise ValueError('Use at least 10 reads per transcript and 100 bp of reference padding')
    generated_names = ['reference.fa', 'reference.fa.fai', 'annotation.gtf', 'annotation.junctions.bed',
                       'reference_regions.tsv', 'reads.fastq.gz', 'inputs.tsv', 'selected_reads.tsv', 'provenance.json']
    if {p.resolve() for p in (args.bam, args.fasta, args.gtf)} & {(args.output/name).resolve() for name in generated_names}:
        raise ValueError('Example outputs must not overwrite source inputs')
    args.output.mkdir(parents=True, exist_ok=True)
    records, gtf_hash = reference_records(args.gtf)
    exon_coordinates = defaultdict(set)
    metadata = {}
    for tid, fields, parsed in records:
        metadata[tid] = parsed
        if fields[2] == 'exon':
            exon_coordinates[tid].add((int(fields[3]), int(fields[4])))
    models = {}
    for tid in TRANSCRIPTS:
        exons = sorted(exon_coordinates[tid])
        if len(exons) < 2:
            raise ValueError(f'Expected a spliced transcript: {tid}')
        parsed = metadata[tid]
        attrs = parsed['attributes']
        models[tid] = dict(transcript_id=tid, chromosome=parsed['seqname'], strand=parsed['strand'],
            gene=attrs.get('gene_name') or attrs['gene_id'], gene_id=attrs.get('ncbi_gene_id', ''),
            start=exons[0][0], end=exons[-1][1], exons=exons,
            intron_chain=tuple(v for left, right in zip(exons, exons[1:]) for v in (left[1], right[0])))
    if len({(m['chromosome'], m['intron_chain']) for m in models.values()}) != len(models):
        raise ValueError('Example transcript models must have distinct complete intron chains')
    regions = {}
    with pysam.FastaFile(str(args.fasta)) as genome:
        for gene in sorted({m['gene'] for m in models.values()}):
            group = [m for m in models.values() if m['gene'] == gene]
            chrom = group[0]['chromosome']
            start = max(0, min(m['start'] for m in group) - 1 - args.flank)
            end = min(genome.get_reference_length(chrom), max(m['end'] for m in group) + args.flank)
            regions[gene] = dict(contig=gene+'_region', original_chromosome=chrom,
                original_start_0based=start, original_end_0based_exclusive=end,
                transcript_ids=';'.join(sorted(m['transcript_id'] for m in group)))
        fasta = args.output/'reference.fa'
        with fasta.open('w') as handle:
            for gene, region in regions.items():
                sequence = genome.fetch(region['original_chromosome'], region['original_start_0based'],
                                        region['original_end_0based_exclusive'])
                handle.write('>'+region['contig']+'\n')
                for i in range(0, len(sequence), 60):
                    handle.write(sequence[i:i+60]+'\n')
        (args.output/'reference.fa.fai').unlink(missing_ok=True)
        pysam.faidx(str(fasta))
    write_tsv(args.output/'reference_regions.tsv', list(regions.values()))
    annotation = args.output/'annotation.gtf'
    written = set()
    with annotation.open('w') as handle:
        handle.write('# RefSeq GRCh38.p14 subset; local contig coordinates; see reference_regions.tsv\n')
        for tid, fields, _ in records:
            region = regions[models[tid]['gene']]
            fields = fields.copy()
            fields[0] = region['contig']
            fields[3] = str(int(fields[3]) - region['original_start_0based'])
            fields[4] = str(int(fields[4]) - region['original_start_0based'])
            if not 1 <= int(fields[3]) <= int(fields[4]) <= region['original_end_0based_exclusive']-region['original_start_0based']:
                raise ValueError('An annotation feature lies outside its cropped reference')
            line = '\t'.join(fields)
            if line not in written:
                handle.write(line+'\n')
                written.add(line)
    converter = ROOT/'trek/gtf2bed.pl'
    bed = subprocess.check_output(['perl', str(converter), str(annotation)], text=True)
    bed_lines = sorted(bed.splitlines(), key=lambda line: (line.split('\t')[0], int(line.split('\t')[1]), line.split('\t')[3]))
    (args.output/'annotation.junctions.bed').write_text('\n'.join(bed_lines)+'\n')
    selected = []
    selection_summary = []
    seen = set()
    with pysam.AlignmentFile(str(args.bam), 'rb') as bam:
        if not bam.has_index():
            raise ValueError('The source BAM must have a BAI or CSI index')
        for tid, model in models.items():
            region = regions[model['gene']]
            eligible, rejected = {}, Counter()
            for read in bam.fetch(model['chromosome'], model['start']-1, model['end']):
                if read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < args.min_mapq:
                    rejected['alignment_filter'] += 1
                    continue
                if _extract_junctions(read.cigartuples, read.reference_start) != model['intron_chain']:
                    rejected['different_intron_chain'] += 1
                    continue
                if read.query_sequence is None or read.query_qualities is None or any(op == 5 for op, _ in read.cigartuples):
                    rejected['unrecoverable_original_fastq'] += 1
                    continue
                if read.reference_start < region['original_start_0based'] or read.reference_end > region['original_end_0based_exclusive']:
                    rejected['outside_reference_crop'] += 1
                    continue
                sequence = read.get_forward_sequence()
                qualities = read.get_forward_qualities()
                if len(sequence) != len(qualities) or any(q < 0 or q > 93 for q in qualities):
                    rejected['invalid_fastq_quality'] += 1
                    continue
                quality = pysam.qualities_to_qualitystring(qualities)
                if read.query_name in eligible:
                    raise ValueError(f'Duplicate eligible primary read name: {read.query_name}')
                eligible[read.query_name] = dict(read_id=read.query_name, transcript_id=tid, gene=model['gene'],
                    original_chromosome=model['chromosome'], original_start_0based=read.reference_start,
                    original_end_0based_exclusive=read.reference_end, original_reverse=read.is_reverse,
                    source_mapq=read.mapping_quality, read_length=len(sequence),
                    sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
                    quality_sha256=hashlib.sha256(quality.encode()).hexdigest(), sequence=sequence, quality=quality)
            order = sorted(eligible, key=lambda name: (hashlib.sha256(f'{args.seed}|{tid}|{name}'.encode()).digest(), name))
            chosen = order[:args.reads_per_transcript]
            if len(chosen) != args.reads_per_transcript:
                raise ValueError(f'{tid}: only {len(chosen)} eligible reads for requested {args.reads_per_transcript}')
            for name in chosen:
                if name in seen:
                    raise ValueError('A read was selected for more than one transcript')
                seen.add(name)
                selected.append(eligible[name])
            selection_summary.append(dict(transcript_id=tid, gene=model['gene'], strand=model['strand'],
                exons=len(model['exons']), eligible_reads=len(eligible), selected_reads=len(chosen),
                excluded_alignment_counts=dict(rejected)))
    selected.sort(key=lambda row: (row['transcript_id'], row['read_id']))
    fastq = args.output/'reads.fastq.gz'
    with fastq.open('wb') as raw:
        with gzip.GzipFile(filename='', fileobj=raw, mode='wb', mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding='ascii', newline='\n') as handle:
                for row in selected:
                    handle.write(f"@{row['read_id']}\n{row['sequence']}\n+\n{row['quality']}\n")
    provenance_rows = [{k:v for k,v in row.items() if k not in ('sequence','quality')} for row in selected]
    write_tsv(args.output/'selected_reads.tsv', provenance_rows)
    (args.output/'inputs.tsv').write_text('fastq\tplatform\nreads.fastq.gz\tdrna\n')
    generated = ['reads.fastq.gz', 'reference.fa', 'reference.fa.fai', 'annotation.gtf',
                 'annotation.junctions.bed', 'reference_regions.tsv', 'selected_reads.tsv', 'inputs.tsv']
    provenance = dict(schema=1, species='Homo sapiens', sample='K562', platform='Oxford Nanopore',
        library_type='direct RNA', run_accession='SRR9304718', bioproject='PRJNA548942',
        public_run_url='https://www.ncbi.nlm.nih.gov/sra/?term=SRR9304718',
        reference_assembly='GCF_000001405.40_GRCh38.p14', source_bam_sha256=sha256(args.bam),
        source_gtf_sha256=gtf_hash, source_fasta_name=args.fasta.name,
        source_bam_mapping=dict(aligner='minimap2 2.30-r1287', preset='splice', splice_strand='b', kmer=14,
            reference='full GRCh38.p14', junction_guide='RefSeq BED12', secondary='no'),
        sequence_policy='Original query sequences and qualities recovered from primary BAM records; reverse-complemented BAM sequences and reversed quality arrays restored to original read orientation. No sequence trimming, synthetic bases, or fabricated quality values.',
        selection=dict(seed=args.seed, min_mapq=args.min_mapq, reads_per_transcript=args.reads_per_transcript,
            rule='Lowest SHA256(seed|transcript_id|read_id) ranks within exact full-intron-chain matches; exclude hard-clipped or quality-missing records and alignments outside the reference crop.',
            sample_is_for_demonstration_not_inference=True, transcripts=selection_summary),
        reference_regions=list(regions.values()),
        coordinates='GTF and output PAS positions are 1-based on cropped contigs. Add original_start_0based from reference_regions.tsv to recover original 1-based genome positions. BED12 is 0-based half-open.',
        build_script_sha256=sha256(Path(__file__)),
        files={name:dict(bytes=(args.output/name).stat().st_size, sha256=sha256(args.output/name)) for name in generated})
    (args.output/'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    print(json.dumps(dict(reads=len(selected), transcripts=selection_summary,
        genomic_reference_bases=sum(r['original_end_0based_exclusive']-r['original_start_0based'] for r in regions.values()),
        data_bytes=sum(provenance['files'][name]['bytes'] for name in generated)), indent=2))


if __name__ == '__main__':
    main()
