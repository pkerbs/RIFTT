# Example

A minimal end-to-end run on three publicly available LL-100 cell lines, each
carrying a well-established fusion:

| sample   | fusion             | ENA run    | ENA sample   |
| -------- | ------------------ | ---------- | ------------ |
| K-562    | `BCR::ABL1`      | ERR3003547 | SAMEA5178166 |
| KASUMI-1 | `RUNX1::RUNX1T1` | ERR3003548 | SAMEA5178167 |
| NB-4     | `PML::RARA`      | ERR3003575 | SAMEA5178194 |

Project [PRJEB30312](https://www.ebi.ac.uk/ena/browser/view/PRJEB30312) (LL-100 panel).

## Run

1. Download the paired FASTQ files for the three runs from ENA, e.g.

   ```
   for r in ERR3003547 ERR3003548 ERR3003575; do
     wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${r:0:6}/00${r: -1}/$r/${r}_1.fastq.gz"
     wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${r:0:6}/00${r: -1}/$r/${r}_2.fastq.gz"
   done
   ```

   and rename them to `<sample>_R1.fastq.gz` / `<sample>_R2.fastq.gz`
   (`K-562`, `KASUMI-1`, `NB-4`).
2. Configure RIFTT (`config/params.conf`) with your reference files.
3. Detection:

   ```
   for s in K-562 KASUMI-1 NB-4; do ./RIFTT_part1.sh "$s" /path/to/fastq; done
   ```
4. Filtering, using the clinical table in this folder:

   ```
   # set CLINTABLE=example/clintable.xlsx in config/params.conf
   ./RIFTT_part2.sh
   ```

Part 2 also needs the GENCODE v38 annotation GTF (`ANNO` in the config).

## Expected result

`BCR::ABL1`, `RUNX1::RUNX1T1` and `PML::RARA` are each called with
`ev_level == 7`. `expected_resultTable.tsv` holds the high-evidence calls.

## Files

| file                         | content                                                                                                        |
| ---------------------------- | -------------------------------------------------------------------------------------------------------------- |
| `clintable.xlsx`           | Part 2 input table for the three samples                                                                       |
| `expected_resultTable.tsv` | the`ev_level >= 6` rows of the reference run                                                                 |
| `part2_intermediates/`     | the Arriba / FusionCatcher / featureCounts / insert-size files Part 2 reads, to test Part 2 without alignment |
