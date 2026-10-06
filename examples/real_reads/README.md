# Real-read example

This example runs FASTQ-to-PAS identification using 300 actual Oxford Nanopore
direct-RNA reads from human K562 run
[SRR9304718](https://www.ncbi.nlm.nih.gov/sra/?term=SRR9304718), BioProject
[PRJNA548942](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA548942).

## Files

| File | Purpose |
| --- | --- |
| `reads.fastq.gz` | Original sequences and quality values recovered from primary BAM records |
| `inputs.tsv` | FASTQ path and library type (`drna`); the path is relative to this file |
| `reference.fa` and `reference.fa.fai` | Two contiguous genomic reference crops, totalling 22,689 bp |
| `annotation.gtf` | Three RefSeq models with coordinates shifted onto the cropped reference |
| `annotation.junctions.bed` | Matching BED12 junction guide |
| `reference_regions.tsv` | Mapping to original GRCh38.p14 coordinates |
| `selected_reads.tsv` | Read names, source alignments, model selection and sequence/quality hashes |
| `provenance.json` | Accessions, reference version, selection rule and input checksums |
| `expected/` | Expected PAS tables, summary, assignment counts/end positions and software versions |
| `check_results.py` | Verification of a completed example run |
| `build_example.py` | Optional reconstruction from the full source BAM and references |

The example is about 0.4 MB. No full-genome FASTA, annotation or BAM is needed
for the quick start.

## Run

After installing TREK with `environment.yml`, run from the repository root:

```bash
trek -g examples/real_reads/annotation.gtf \
     -f examples/real_reads/reference.fa \
     -i examples/real_reads/inputs.tsv \
     -o example_results -p example \
     -t 2 --serial --min-overlap 0.5

python examples/real_reads/check_results.py --output example_results
```

Use a fresh output directory for a new end-to-end run. An existing output
directory/prefix reuses its read-assignment pickle only when the accompanying
`.read_assignments.filters.json` confirms mapping completion, SA-tag filtering,
matching library/assignment settings and matching input/reference paths.
The input list explicitly uses the direct-RNA profile (`splice -u f -k14`).
Examples run with an older TREK version should use a new output directory or prefix.
The expected outputs correspond to this supplied input list. Omitting the type
uses the generic `splice -u b -k15` profile and can produce different assignments
and PAS support counts.

## Expected results

| Transcript | Gene | Strand | Exons | Assigned reads | Retained PASs |
| --- | --- | --- | ---: | ---: | ---: |
| NM_001289746.2 | GAPDH | + | 8 | 100 | 1 |
| NM_002046.7 | GAPDH | + | 9 | 100 | 1 |
| NM_181697.3 | PRDX1 | - | 6 | 100 | 2 |

The GAPDH models have different internal splice chains and share the reported
PAS at `GAPDH_region:5855`. PRDX1 has sites at `PRDX1_region:2003` and
`PRDX1_region:2177`, with support 72 and 28. The output therefore has four
transcript-PAS records at three distinct genomic positions. No site is removed
by internal-priming filtering here, so that output contains only its header.

GAPDH PAS support counts are 60 and 76, respectively, despite 100 assigned
reads per model. Support from components rejected by the peak filters is not
added to the surviving PAS. Exact reference outputs are in `expected/`.

## Coordinates and selection

Reference crops preserve the original genomic sequence and intron lengths.
Transcript structures and RefSeq IDs are retained; only genomic coordinates and
contig names are translated to keep the reference small.

| Example contig | Original interval (0-based, half-open) |
| --- | --- |
| GAPDH_region | NC_000012.12:6532516-6540371 |
| PRDX1_region | NC_000001.11:45509050-45523884 |

Add `original_start_0based` from `reference_regions.tsv` to a local 1-based PAS
position to recover its original 1-based genomic coordinate. For example,
`GAPDH_region:5855` corresponds to `NC_000012.12:6538371`.

The source BAM was aligned to the full RefSeq GRCh38.p14 genome. Read selection
requires a primary alignment with MAPQ >=20, an exact complete intron-chain
match, recoverable sequence/qualities and an alignment contained in the crop.
Hard-clipped reads are excluded. Within each model, 100 reads were selected by
deterministic SHA-256 ranking using seed `20261004`.

Reverse-aligned BAM records were restored to their original sequencing
orientation, including quality order. No bases or quality values were invented
or trimmed. The models were chosen to show shared PASs across splice structures
and multiple PASs on a negative-strand transcript. The equal read counts and
reduced reference are for demonstrating execution, not estimating full-sample
expression or PAS usage.

## Rebuild from source data

Rebuilding is optional and requires the original SRR9304718 BAM with its BAI/CSI
index and the full references used for mapping:

```bash
python examples/real_reads/build_example.py \
  --bam /path/to/SRR9304718/alignments.bam \
  --fasta /path/to/GCF_000001405.40_GRCh38.p14_genomic.fna \
  --gtf /path/to/GCF_000001405.40_GRCh38.p14_genomic.gtf
```

The builder regenerates the inputs and BED12 guide with the same model list and
seed. It does not overwrite `expected/`.

Reference assembly:
[GCF_000001405.40_GRCh38.p14](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/).
Source BAM/GTF checksums, alignment settings and per-read evidence are retained
in the provenance files. Exact output comparisons use the tested environment
recorded in `expected/versions.json`.
