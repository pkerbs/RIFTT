# RIFTT

**R**obust **I**dentification of **F**usion genes in **T**umor **T**ranscriptomes.

RIFTT calls fusion genes from RNA-seq data with Arriba and FusionCatcher and applies knowledge-based filtering and three scoring metrics (Promiscuity Score, Fusion Transcript Score, Robustness Score) to prioritise robust candidates. It is bundled in a Singularity container and runs in two parts: detection per sample (Part 1) and cohort-level filtering (Part 2).

## Requirements

- Singularity or Apptainer, version 3.6 or newer
- at least 40 GB RAM
- Reference files:

  - genome FASTA and gene annotation GTF (GENCODE, [https://www.gencodegenes.org](https://www.gencodegenes.org))
  - STAR index built **without** `--sjdbGTFfile` / `--sjdbOverhang`
  - FusionCatcher genomic database
    ([https://sourceforge.net/projects/fusioncatcher/files/data](https://sourceforge.net/projects/fusioncatcher/files/data))

## Build the container

```
git clone https://github.com/pkerbs/RIFTT.git
cd RIFTT
sudo singularity build RIFTT.sif RIFTT.def       # or: singularity build --fakeroot
```

`RIFTT.sif` is expected next to the helper scripts (or set `RIFTT_SIF`).

## Configure

```
cp config/params.conf.example config/params.conf
```

Edit `config/params.conf` (paths, genome build, threads, ...). Both helper scripts read this file. All options are described in `config/params.conf.example`.

## Run

**Part 1 - detection (once per sample):**

```
./RIFTT_part1.sh <sample_name> <fastq_folder>
```

FASTQ files are matched as `<sample_name>*R1.fastq.gz` / `*R2.fastq.gz`.

**Part 2 - filtering (once, over the whole cohort):**

```
./RIFTT_part2.sh
```

## Clinical table (Part 2 input)

An `.xlsx` file with one row per sample and four columns:

| column        | required | content                                                                                                             |
| ------------- | -------- | ------------------------------------------------------------------------------------------------------------------- |
| `cohort`    | yes      | cohort / group label; PS, maxPS and RS are computed over all samples of the run, and the plots are drawn per cohort |
| `sample`    | yes      | sample name, identical to the one used in Part 1                                                                    |
| `Karyotype` | no       | ISCN karyotype string, e.g.`46,XX,t(6;9)(p23;q34)[12]`; leave empty if unknown                                    |
| `otherCyto` | no       | fusions known for this sample from other diagnostics; leave empty if none                                           |

`Karyotype` and `otherCyto` are annotation only - they fill the `karyo` / `mol`
columns of the output and the oncoprint, and do **not** affect filtering or
`ev_level`. A clinical table with only `cohort` and `sample` filled in is valid.

`otherCyto` holds comma-separated entries, each one of:

| entry            | meaning                                     | sets`mol` |
| ---------------- | ------------------------------------------- | ----------- |
| `GENE1::GENE2` | fusion confirmed by other diagnostics       | `Y`       |
| `GENE#`        | one partner screened / suspected            | `S`       |
| `GENE1!GENE2`  | fusion explicitly tested and found negative | `N`       |

Example values: `PML::RARA` · `RUNX1::RUNX1T1,CDC2L1::SLC35E2` · `DEK!NUP214`

## Output

Part 2 writes to `<OUTPUTFOLDER>/filter_results/run_<timestamp>/`:

| file                                                | content                                             |
| --------------------------------------------------- | --------------------------------------------------- |
| `resultTable.xlsx`                                | all evaluated fusion calls with their metrics       |
| `filterrun.RData`                                 | the R workspace of the run                          |
| `plots/Violinplot_PS.png`, `Violinplot_FTS.png` | PS / FTS distributions, known vs unknown            |
| `plots/TPM-FTS_3D_plot.html`                      | interactive TPM vs FTS scatter                      |
| `plots/Oncoprint_Karyo_MDx_RNAseq.png`            | oncoprint of high-evidence fusions vs clinical data |
| `plots/*_circos.png`                              | circos plot of robust fusions per cohort            |

### `resultTable.xlsx` columns

| column                                                            | meaning                                                                                          |
| ----------------------------------------------------------------- | ------------------------------------------------------------------------------------------------ |
| `cohort`, `sample`                                            | from the clinical table                                                                          |
| `caller`                                                        | `AR` (Arriba) or `FC` (FusionCatcher)                                                        |
| `gene1`, `gene2`, `label`                                   | 5' and 3' partner gene and the`gene1::gene2` label                                             |
| `reciprocal`                                                    | `TRUE` if this is the reciprocal orientation of a call also reported the other way round       |
| `break5prime`, `break3prime`                                  | breakpoint coordinates                                                                           |
| `cov`                                                           | breakpoint-spanning read count                                                                   |
| `mitelman_rec`                                                  | number of MitelmanDB records for this fusion                                                     |
| `PS`                                                            | Promiscuity Score (mean number of distinct partners of the two genes)                            |
| `tpm5prime`, `tpm3prime`, `tpmfusion`                       | TPM of the 5' partner, 3' partner and the fusion transcript                                      |
| `FTS5`, `FTS3`, `FTS`                                       | Fusion Transcript Score of the 5' side, 3' side, and their mean                                  |
| `RS`                                                            | Robustness Score (fraction of the samples calling this fusion in which FTS passes)               |
| `known`                                                         | `known` if the fusion is a validated ChimerDB entry (>= 2 PMIDs, PCR/Sanger), else `unknown` |
| `passCallFilt`, `passBL`, `passPS`, `passFTS`, `passRS` | whether each criterion is met                                                                    |
| `callerOverlap`                                                 | `TRUE` if called by both Arriba and FusionCatcher                                              |
| `ev_level`                                                      | evidence level 0-7 (see below)                                                                   |
| `karyo`, `mol`                                                | agreement with the`Karyotype` / `otherCyto` column (`Y` / `S` / `N`)                   |

### Evidence level

`ev_level` counts how many of these seven criteria a call meets:

1. `passCallFilt` - passes the caller's built-in filters
2. `passBL` - not on the healthy-tissue blacklist
3. `passPS` - Promiscuity Score <= the highest PS among known fusions in the run (`maxPS`)
4. `passFTS` - FTS5, FTS3 >= 0.025; FTS5, FTS3 < 1; FTS >= 0.05
5. `passRS` - Robustness Score >= 0.5
6. `callerOverlap` - called by both callers
7. `known` - validated ChimerDB entry

`ev_level >= 6` marks high-confidence candidates.

## Example

Example on three public LL-100 cell lines: see `example/`.
