# Stage 1 — Build the non-redundant TE consensus library

This stage produces the two inputs TEMP2 needs:
- `LTRRT.consensus2.fa` — non-redundant TE consensus library
- `ref.fa.out2.bed6` — reference TE annotation coordinates (BED6)

It uses **synLTR** (`module2`): https://github.com/cwb14/synLTR @ `990df32`.
Clone synLTR at that commit and run from `synLTR/module2/`. Use the conda
environment in `data/synLTR_ltrharvest_env.yml`.

> Note on script names: in the version of synLTR used here, the consensus step
> is `consensus_library.sh` (it was previously `consensus_ltrrt_library2.sh`).
> `ltrharvest.py` and `mask_ltr.py` are unchanged.

The procedure builds LTR-RT annotations iteratively, masking each round's hits
before the next so that nested elements are resolved (un-nested → single → double
→ triple nest), then merges everything into one consensus library.

## Genome

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/010/389/155/GCF_010389155.1_C_canadensis_v1/GCF_010389155.1_C_canadensis_v1_genomic.fna.gz
# renamed to ref.fa below
```

## Iterative LTR-RT discovery

```bash
# Un-nested LTR-RTs
python ltrharvest.py --genome ref.fa --proteins ref.pep --threads 120 --out-prefix ref_ltr \
  --scn-min-ltr-len 100 --scn-min-ret-len 800 --scn-max-ret-len 15000 \
  --scn-min-int-len 500 --scn-max-int-len 12000

# Single-nest LTR-RTs
python mask_ltr.py --features-fasta ref_ltr.ltrharvest.full_length.dedup.fa.rexdb-plant.cls.lib.fa \
  --genome ref.fa --feature-character N --far-character V --distance 15000 > ref_r1.fa
python ltrharvest.py --require-run-chars N --genome ref_r1.fa --proteins ref.pep --threads 200 \
  --out-prefix ref_r2 --scn-min-ltr-len 100 --scn-min-ret-len 1000 --scn-max-ret-len 30000 \
  --scn-min-int-len 500 --scn-max-int-len 28000 \
  --ltrharvest-args '-mindistltr 100 -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -similar 70 -vic 30 -seed 15 -seqids yes -xdrop 10 -maxdistltr 30000' \
  --ltrfinder-args  '-w 2 -C -D 30000 -d 100 -L 7000 -l 100 -p 20 -M 0.00 -S 0.0'

# Double-nest LTR-RTs
python mask_ltr.py --features-fasta ref_r2.ltrharvest.full_length.dedup.fa.rexdb-plant.cls.lib.fa \
  --genome ref.fa --feature-character R --far-character V --distance 15000 > ref_r2.fa
python ltrharvest.py --require-run-chars N,R --genome ref_r2.fa --proteins ref.pep --threads 100 \
  --out-prefix ref_r3 --scn-min-ltr-len 100 --scn-min-ret-len 1000 --scn-max-ret-len 45000 \
  --scn-min-int-len 500 --scn-max-int-len 43000 \
  --ltrharvest-args '-mindistltr 100 -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -similar 70 -vic 30 -seed 15 -seqids yes -xdrop 10 -maxdistltr 45000' \
  --ltrfinder-args  '-w 2 -C -D 45000 -d 100 -L 7000 -l 100 -p 20 -M 0.00 -S 0.0'

# Triple-nest LTR-RTs
python mask_ltr.py --features-fasta ref_r3.ltrharvest.full_length.dedup.fa.rexdb-plant.cls.lib.fa \
  --genome ref.fa --feature-character D --far-character V --distance 15000 > ref_r3.fa
python ltrharvest.py --require-run-chars N,R,D --genome ref_r3.fa --proteins ref.pep --threads 100 \
  --out-prefix ref_r4 --scn-min-ltr-len 100 --scn-min-ret-len 1000 --scn-max-ret-len 50000 \
  --scn-min-int-len 500 --scn-max-int-len 50000 \
  --ltrharvest-args '-mindistltr 100 -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -similar 70 -vic 30 -seed 15 -seqids yes -xdrop 10 -maxdistltr 50000' \
  --ltrfinder-args  '-w 2 -C -D 50000 -d 100 -L 7000 -l 100 -p 20 -M 0.00 -S 0.0'
```

## Whole-genome TE annotation with EDTA

In parallel with the iterative LTR-RT discovery above, the genome was annotated
with **EDTA** (Ou et al. 2019) to capture intact elements across all TE classes.
The `*.EDTA.intact.fa` output is merged with the LTR-RT rounds below.

```bash
perl EDTA.pl --genome ref.fa --cds ref.cds --curatedlib TE.lib \
  --overwrite 1 --sensitive 1 --anno 1 --threads 40
# → ref.fa.mod.EDTA.intact.fa  (intact TEs across classes)
```

Inputs:
- `ref.cds` — coding sequences for the assembly, passed to EDTA to purge
  gene-related false positives (`--cds`).
- `TE.lib` — curated TE library supplied to EDTA (`--curatedlib`). This study
  used **RepBase** (Bao et al. 2015). RepBase is distributed by GIRI under a
  separate license/subscription, so it is not redistributed here; obtain it from
  https://www.girinst.org/repbase/ .

EDTA: https://github.com/oushujun/EDTA — Ou et al. (2019) *Genome Biology* 20:275.

## Merge rounds into one library

```bash
# Remove masked bases from nested rounds
cat ref*.ltrharvest.full_length.dedup.fa.rexdb-plant.cls.lib.fa \
  | sed '/^>/! { s/[^ATCGatcg]//g; /^$/d }' > ref_all.fa
# Clean header characters that break consensus building
sed '/^>/ s/#.*//' ref_all.fa > ref_all2.fa
cat ref_all2.fa | sed 's/:/_/g' > ref_all3.fa
# Add EDTA intact elements
cat ref.fa.mod.EDTA.intact.fa ref_all3.fa | sed '/^>/ s/#.*//' | sed 's/:/_/g' > ref_all4.fa
```

## Build the consensus library

```bash
bash ./consensus_library.sh -i ref_all4.fa -o TE_consensus_out4 -p LTRRT -t 232 \
  --min-seq-id 0.80 -c 0.80 --cov-mode 0 --cluster-mode 0 --mafft-mode auto --fancy
# → TE_consensus_out4/LTRRT.consensus.fa   (== LTRRT.consensus2.fa used downstream)
```

## Reference TE annotation (BED6 for TEMP2 `-t`)

```bash
# Annotate the genome with the consensus library
RepeatMasker -pa 70 -lib TE_consensus_out4/LTRRT.consensus.fa -no_is ref.fa

# Convert RepeatMasker .out to BED6 (drop headers/blank lines and simple repeats)
awk 'BEGIN{OFS="\t"}
     NR<=2 || $0 ~ /^[[:space:]]*$/ {next}
     $11=="Simple_repeat" {next}
     {
       chr=$5; start=$6-1; end=$7; name=$10; score=0;
       strand=($9=="C" ? "-" : $9);
       print chr, start, end, name, score, strand
     }' ref.fa.out > ref.fa.out2.bed6
```

The two products (`LTRRT.consensus2.fa`, `ref.fa.out2.bed6`) are the `-R` and
`-t` inputs to TEMP2 in stage 2.
