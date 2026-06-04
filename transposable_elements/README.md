# Transposable element (TE) analyses

Scripts and workflow for the transposable-element analyses in the marestail
(*Erigeron canadensis*) herbicide-resistance manuscript.

The TE analysis has four stages. Two of them are run with standalone tools that
are maintained in their own repositories (and cited in the manuscript); this
directory does **not** copy those tools, it documents exactly how they were run
and pins the versions used. The genuinely project-specific glue — the
annotation step and the figure/statistics generator — lives here in full.

```
                                                    tool / script              this repo?
1. Build a non-redundant TE consensus library   →  synLTR + EDTA               referenced
2. Detect TE insertions across 345 samples       →  TEMP2 (fork)               referenced + driver here
3. Annotate insertions with TE family            →  add_transposon_family.py   YES
4. Figures + statistics                           →  te_figures.py             YES
```

## External tools (cite + install these)

| Tool | Used for | Version pinned for this study | License |
|------|----------|-------------------------------|---------|
| **synLTR** | Stage 1 — iterative LTR-RT discovery | `cwb14/synLTR` @ `990df32` | see repo |
| **EDTA** | Stage 1 — whole-genome TE annotation (merged with LTR-RTs) | Ou et al. 2019 | GPL-3.0 |
| **RepBase** | Stage 1 — curated TE library passed to EDTA (`--curatedlib`) | GIRI (licensed) | proprietary |
| **TEMP2** | Stage 2 — insertion detection | `cwb14/TEMP2` @ `45cd113` (fork of `weng-lab/TEMP2`) | GPL-3.0 |

- synLTR: https://github.com/cwb14/synLTR  (commit `990df32acbb55581ccd1978438ebec787a123d96`)
- EDTA: https://github.com/oushujun/EDTA — Ou et al. (2019) *Genome Biology* 20:275.
- RepBase: https://www.girinst.org/repbase/ — Bao et al. (2015) *Mobile DNA* 6:11 (licensed; not redistributed here).
- TEMP2 (this study's fork, with minor bug-fixes): https://github.com/cwb14/TEMP2 (commit `45cd113892db0c27cc450a0a6401444007949d58`)
- TEMP2 upstream / primary citation: Yu et al. (2021) *Nucleic Acids Research* — "A benchmark and an algorithm for detecting germline transposon insertions and measuring de novo transposon insertion frequencies." https://github.com/weng-lab/TEMP2

> The TEMP2 fork differs from upstream only in small bug-fixes to
> `TEMP2_insertion2.sh` and `pickUniqPairFastq.sh`; pin to the commit above to
> reproduce exactly.

## Layout

```
transposable_elements/
├── README.md                     ← you are here
├── environment.yml               ← conda env for stages 3–4 (annotation + figures)
├── 01_build_te_library.md        ← exact synLTR commands (stage 1)
├── 02_run_temp2/
│   ├── run_TEMP2.slurm           ← SLURM array driver (one task per sample)
│   └── prefixes.txt              ← the 345 sample prefixes
├── 03_annotate_families/
│   └── add_transposon_family.py  ← adds TE family to each insertion
├── 04_figures/
│   └── te_figures.py             ← all manuscript TE figures + stats
└── data/
    ├── README.md                 ← how to obtain the large input/derived files
    └── synLTR_ltrharvest_env.yml ← exact conda env used for stage 1
```

## Reproducing the analysis

### Stage 1 — TE consensus library
See [`01_build_te_library.md`](01_build_te_library.md). Produces the consensus TE
library (`LTRRT.consensus2.fa`) and reference TE annotations
(`ref.fa.out2.bed6`) that feed TEMP2. Run with the synLTR environment
(`data/synLTR_ltrharvest_env.yml`).

### Stage 2 — TEMP2 insertion detection
`02_run_temp2/run_TEMP2.slurm` is a SLURM array job that runs `TEMP2 insertion2`
once per sample listed in `prefixes.txt` (345 tasks, 10 concurrent). Before
submitting, edit the input-path variables at the top of the script and the
`conda activate TEMP2` line for your environment. Each sample produces
`<prefix>_TEMP2/<prefix>.insertion.bed`.

```bash
sbatch 02_run_temp2/run_TEMP2.slurm
```

### Stage 3 — Annotate TE families
Adds the TE family (from the consensus library headers) to each insertion call:

```bash
for d in *_TEMP2; do
  prefix="${d%_TEMP2}"
  python 03_annotate_families/add_transposon_family.py \
      "$d/$prefix.insertion.bed" family_list.tsv > "$d/$prefix.insertion.fam.bed"
done
```

`family_list.tsv` is derived from the consensus-library FASTA headers
(`>chr:start-end#LTR/Ty1/SIRE` style); see `data/README.md`.

### Stage 4 — Figures and statistics
`04_figures/te_figures.py` builds the five-page multi-panel figure PDF and the
summary statistics from the per-sample annotated BED files plus a master sample
table. Recommended command (marestail dataset):

```bash
python 04_figures/te_figures.py \
    --master      MASTER_MARESTAIL_MERGED.tsv \
    --te-pattern  '{sample}_TEMP2/{sample}.insertion.fam.bed' \
    --output      te_figures.pdf \
    --awk         '$5 >= 0.1 && $8 >= 3 && $7 == "1p1"'
```

The `--awk` filter keeps high-confidence calls: insertion frequency ≥ 10 %
(`$5`), ≥ 3 supporting reads (`$8`), and split-read evidence at both TE ends
(`$7 == "1p1"`). Run `python 04_figures/te_figures.py -h` for the full set of
optional inputs (`--fai`, `--gff`, `--crm`, `--ltr-age`, `--mutation-rate`,
`--repeatmasker-tbl`, `--sample-col`).

## Environments
- Stages 3–4: `environment.yml` (this directory).
- Stage 1: `data/synLTR_ltrharvest_env.yml` (exact env used; large, fully pinned).
