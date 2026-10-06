# TREK

**Transcript-resolved polyadenylation sites from long-read RNA-seq.**

TREK identifies polyadenylation sites (PASs, transcription end sites) from
Oxford Nanopore and PacBio RNA-seq reads and reports every site on the annotated
transcript it was observed on. Reads are assigned to transcript models by their
complete splice structure before their 3′ ends are clustered, so a PAS is never
inferred for an isoform it was not seen on.

Three commands cover the workflow:

| Command | Scope | What it does |
| --- | --- | --- |
| `trek` | one sample | align reads, assign them to transcripts, cluster read 3′ ends into PASs, remove internal-priming artefacts |
| `unify_apa_sites` | one species | merge per-sample PAS tables into a non-redundant catalogue with stable cluster IDs |
| `annotate_apa` | one species | add sequence context, lineage-specific PAS motifs and PAS-L1–L3 evidence levels |

```
FASTQ + genome + annotation
        │  trek  (per sample)
        ▼
<sample>.apa_sites.txt ──┐
<sample>.apa_sites.txt ──┤  unify_apa_sites  (per species)
<sample>.apa_sites.txt ──┘
        │
        ▼
<species>_unified_apa_sites.txt
        │  annotate_apa
        ▼
<species>_unified_apa.anno.txt   (motif, PAS-L1/L2/L3, ±50 nt context)
```

## Installation

```bash
git clone https://github.com/TFC-Lab-CUHK/TREK.git
cd TREK
mamba env create -f environment.yml   # Python, minimap2 and pinned dependencies
mamba activate TREK
pip install -e .
```

With an existing Python ≥ 3.12 environment, `pip install -e .` is enough;
install [minimap2](https://github.com/lh3/minimap2) separately and keep it on
`PATH`. Dependencies are listed in [pyproject.toml](pyproject.toml); the pinned
versions in [environment.yml](environment.yml) are the ones the bundled example
was validated with.

## Quick start

A complete run on 300 real direct-RNA reads takes under a minute and needs no
genome download:

```bash
trek -g examples/real_reads/annotation.gtf \
     -f examples/real_reads/reference.fa \
     -i examples/real_reads/inputs.tsv \
     -o example_results -p example -t 2 --serial

python examples/real_reads/check_results.py --output example_results
# PASS: 300 input reads, 3 transcript models; assignments and PAS outputs match.
```

Then merge and annotate two small synthetic samples:

```bash
unify_apa_sites examples/unification/Demo_species --output-dir example_unified

annotate_apa --unified-file example_unified/Demo_species_unified_apa_sites.txt \
             --fasta examples/unification/reference.fa \
             --motif-model LAND_PLANT \
             --output example_unified/Demo_species_unified_apa.anno.txt
```

Expected results are described in [examples/real_reads](examples/real_reads/README.md)
and [examples/unification](examples/unification/README.md).

## `trek` — per-sample PAS calling

```bash
trek -g annotation.gtf -f genome.fa -i inputs.tsv -o results -p sample -t 16
```

**Inputs**

- `-g` annotation in GTF, GFF or GFF3. Multi-exon reads are matched to these
  models by their exact intron chain, so the annotation defines which
  transcripts can receive reads.
- `-f` uncompressed genome FASTA with the same sequence names as the annotation.
- `-i` input list: one FASTQ per row (gzip accepted), optionally followed by the
  library type. All FASTQs in one list are pooled into one sample.

```text
fastq	platform
/data/run1.fastq.gz	drna
/data/run2.fastq.gz	pcr-cdna
/data/run3.fastq.gz
```

Rows are tab-separated (whitespace-separated also works; quote paths containing
spaces). Relative paths resolve from the list's directory. An empty or `NA`
platform uses the generic profile; unknown labels are errors.

Accepted platform labels are `drna` (ONT direct RNA), `pcr-cdna` and `dcdna`
(ONT cDNA) and `pacbio` (Iso-Seq / HiFi / CCS); each selects a matching minimap2
profile. `--splice-strand f` can be given when cDNA or PacBio reads are known to
be in transcript orientation.

**Options you are most likely to change**

| Option | Default | Meaning |
| --- | --- | --- |
| `-t / --threads` | 8 | minimap2 threads |
| `-j / --jobs` | all CPUs | parallel workers for PAS calling; `--serial` forces one |
| `--min-mapq` | 1 | minimum mapping quality of a primary alignment |
| `--min-overlap` | 0.5 | fraction of an annotated single-exon transcript a read must cover |
| `--min-reads` | 10 | reads a transcript needs before its 3′ ends are analysed |
| `--min-cluster-size` | 10 | reads a candidate peak needs to be retained |
| `--max-clusters` | 5 | largest number of mixture components tried |
| `--min-distance` | 50 | minimum spacing (nt) between retained PASs |
| `--min-dominance` | 0.3 | a peak needs at least this fraction of the largest peak's reads |
| `--min-sharpness` | 0.5 | local concentration a peak must show (`1 − IQR / min-distance`) |
| `--priming-window` / `--priming-a-threshold` | 20 / 0.5 | downstream window and A-fraction for the internal-priming filter; `--no-filter-priming` disables it |

**Outputs** (`<output>/<prefix>.*`)

| File | Content |
| --- | --- |
| `apa_sites.txt` | one row per transcript–PAS: coordinate (1-based), supporting reads, within-transcript abundance |
| `internal_priming_removed.txt` | sites removed by the internal-priming filter, with their downstream A-content |
| `summary.txt` | transcripts analysed and the distribution of PASs per transcript |
| `read_assignments.pkl` / `.filters.json` | cached read 3′ ends per transcript and the settings that produced them |

## Examples

- [examples/real_reads](examples/real_reads/README.md) — provenance and expected results of the real-read example
- [examples/unification](examples/unification/README.md) — unification and annotation walkthrough

## Citation

Citation to be added.

## License

MIT
