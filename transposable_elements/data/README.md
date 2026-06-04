# Data files

Large inputs and derived files are **not committed** to keep the repository
lean. This file documents what they are and how to obtain or regenerate them.

## Reference genome
*Erigeron canadensis* assembly **GCF_010389155.1** (C_canadensis_v1), from NCBI:

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/010/389/155/GCF_010389155.1_C_canadensis_v1/GCF_010389155.1_C_canadensis_v1_genomic.fna.gz
```

Used throughout as `ref.fa` (decompressed). The BWA index (`bwa_index/ref.fa`)
referenced by the SLURM driver is built with `bwa index`.

## Derived files (produced by stage 1; regenerate via ../01_build_te_library.md)

| File | What it is | TEMP2 flag |
|------|-----------|------------|
| `LTRRT.consensus2.fa` | Non-redundant TE consensus library | `-R` |
| `ref.fa.out2.bed6`    | Reference TE annotation coordinates (BED6) | `-t` |

These were generated with synLTR (`cwb14/synLTR` @ `990df32`); see
`../01_build_te_library.md` for the exact commands. They total ~35 MB
gzipped — if you prefer to distribute them directly rather than regenerate,
deposit them in a data archive (e.g. Zenodo) and link the DOI here:

> Consensus library + annotations: <ADD ZENODO DOI / RELEASE LINK>

## EDTA inputs (stage 1)
- `ref.cds` — coding sequences for the *E. canadensis* assembly, passed to EDTA
  via `--cds` to reduce gene-derived false positives.
- `TE.lib` — curated TE library passed to EDTA via `--curatedlib`. This study
  used **RepBase** (https://www.girinst.org/repbase/), which is licensed by GIRI
  and therefore not redistributed here.

## Sequencing reads
345 whole-genome short-read samples; accessions listed in the manuscript.
Sample prefixes are in `../02_run_temp2/prefixes.txt`. The SLURM driver expects
paired FASTQs named `<prefix>_1.fastq.gz`/`<prefix>_2.fastq.gz` or
`<prefix>_paired_1.fq.gz`/`<prefix>_paired_2.fq.gz`.

## Annotation / analysis tables
- `family_list.tsv` — maps consensus-library elements (by genomic coordinate) to
  TE family. Derived from the consensus-library FASTA headers
  (`>chr:start-end#LTR/Ty1/SIRE`). Used by `../03_annotate_families/add_transposon_family.py`.
- `MASTER_MARESTAIL_MERGED.tsv` — master sample table (sample IDs, bioclim
  variables, glyphosate phenotypes, etc.) consumed by
  `../04_figures/te_figures.py`. Document its columns alongside the manuscript
  supplementary data.

## synLTR environment
`synLTR_ltrharvest_env.yml` is the exact, fully pinned conda environment used for
stage 1 (genometools/ltrharvest, ltr_finder, mmseqs2, mafft, trimal, etc.).
