# PAS unification example

These two small **synthetic** tables demonstrate the single-sample TREK output
format. The transcript IDs, genes and PAS positions are illustrative, not
biological observations. Both files use the same chromosome naming and coordinate
system. The minus-strand example also uses 1-based genomic coordinates.

From the repository root, run the standalone script:

```bash
python trek/unify_apa_sites.py examples/unification/Demo_species \
    --output-dir example_unified
```

After installing TREK, the equivalent command is:

```bash
unify_apa_sites examples/unification/Demo_species --output-dir example_unified
```

Expected results:

- Three retained genomic clusters, with representative positions 1000, 1300 and
  2000; their inclusive intervals are 1000–1004, 1300–1303 and 1995–2000.
- Ten detailed records retaining their original positions and read counts.
- Three represented transcripts, T1–T3, and two input samples per retained cluster.
- Two filtered-out clusters: G3 at 3000 and G2 at 2300. Each has only one
  transcript–sample observation and fewer than 50 reads.
- T3 in root has relative abundance 1.0 at its retained cluster after filtering.

The detailed table preserves individual transcript–sample–position records;
`unified_ID` and `mode_site_position` link nearby positions to their shared cluster.
The summary table contains one row per genomic cluster within a gene and strand.

## Annotate the unified clusters

The bundled `reference.fa` is a **4,000-bp synthetic sequence**, not the human
chromosome indicated by the demonstration sequence identifier. It places primary
motifs upstream of the example PASs in transcript direction, including the
negative-strand cluster. It is for testing this example only.

```bash
python trek/annotate_apa.py \
    --unified-file example_unified/Demo_species_unified_apa_sites.txt \
    --fasta examples/unification/reference.fa \
    --motif-model LAND_PLANT \
    --output example_unified/Demo_species_unified_apa.anno.txt
```

Expected annotation: three rows, each with `AATAAA` at position −20,
`motif_model=land_plant`, `n_samples=2`, and `pa_level=PAS-L1`. Each displayed
sequence is 101 bp long and oriented in transcript direction. The explicit model
is used because `Demo_species` is synthetic and has no entry in the species table.
