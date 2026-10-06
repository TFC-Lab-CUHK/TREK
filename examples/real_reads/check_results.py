#!/usr/bin/env python3
"""Check bundled inputs and outputs from a completed real-read example run."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import pickle

EXAMPLE = Path(__file__).resolve().parent


def read_table(path):
    with path.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        rows = list(reader)
        return reader.fieldnames, sorted(tuple(row[field] for field in reader.fieldnames) for row in rows)


def check(output, prefix='example'):
    provenance = json.loads((EXAMPLE/'provenance.json').read_text())
    for name, expected in provenance['files'].items():
        if hashlib.sha256((EXAMPLE/name).read_bytes()).hexdigest() != expected['sha256']:
            raise ValueError(f'Example input changed: {name}')
    with gzip.open(EXAMPLE/'reads.fastq.gz', 'rt') as handle:
        n_reads = 0
        while header := handle.readline():
            sequence, plus, quality = handle.readline(), handle.readline(), handle.readline()
            if not header.startswith('@') or plus.strip() != '+' or len(sequence.rstrip()) != len(quality.rstrip()):
                raise ValueError('Malformed bundled FASTQ')
            n_reads += 1
    if n_reads != sum(row['selected_reads'] for row in provenance['selection']['transcripts']):
        raise ValueError('Unexpected read count')
    for suffix in ['apa_sites.txt', 'internal_priming_removed.txt']:
        if read_table(output/f'{prefix}.{suffix}') != read_table(EXAMPLE/'expected'/f'example.{suffix}'):
            raise ValueError(f'{suffix} differs from reference output; check expected/versions.json and the documented command')
    if (output/f'{prefix}.summary.txt').read_text() != (EXAMPLE/'expected/example.summary.txt').read_text():
        raise ValueError('Summary differs from reference output')
    with (EXAMPLE/'expected/assignment_counts.tsv').open() as handle:
        expected = {row['transcript_id']: int(row['assigned_reads']) for row in csv.DictReader(handle, delimiter='\t')}
    # Assignment cache written by TREK during this example run.
    with (output/f'{prefix}.read_assignments.pkl').open('rb') as handle:
        assignments = pickle.load(handle)
    if {tid: len(ends) for tid, ends in assignments.items()} != expected:
        raise ValueError('Transcript assignment counts differ from reference output')
    with (EXAMPLE/'expected/assignment_ends.tsv').open() as handle:
        expected_ends = {row['transcript_id']: [int(pos) for pos in row['positions'].split(',')]
                         for row in csv.DictReader(handle, delimiter='\t')}
    if {tid: sorted(int(pos) for pos in ends) for tid, ends in assignments.items()} != expected_ends:
        raise ValueError('Assigned read-end coordinates differ from reference output')
    print(f'PASS: {n_reads} input reads, {len(expected)} transcript models; assignments and PAS outputs match.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--prefix', default='example')
    args = parser.parse_args()
    check(args.output, args.prefix)


if __name__ == '__main__':
    main()
