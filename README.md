# Introduction to RNA-seq analysis

A 25-minute introduction to bulk RNA-seq, built on one real experiment from raw reads to biology,
plus the pipeline that produces every figure in it. One command runs everything hands-off:

```
FASTQ → quality control → trimming → alignment → count matrix          (upstream, command line)
      → exploration → differential expression → GO / pathways / diseases   (downstream, R)
      → presentation slides
```

**The experiment:** airway smooth muscle cells from 4 donors, each untreated and treated with
dexamethasone (Himes *et al.* 2014, GEO GSE52778; the Bioconductor `airway` dataset). The pipeline
streams the first 2 million read pairs of each of the 8 samples from the European Nucleotide Archive,
so a full run takes about 15 minutes instead of hours.

**The slides:** online at <https://mortezasabri.github.io/rnaseq-workshop/>, or download
[`docs/rnaseq_intro.html`](docs/rnaseq_intro.html) and open it in a web browser (one file, works
offline). Press **F** for full screen.

| Folder | What is in it |
|---|---|
| `run_all.sh` | Runs everything, or one part: `upstream`, `downstream`, `slides` |
| `config/` | `config.sh` (all settings), `samples.tsv` (the 8 airway samples), `samples.template.tsv` (for your own data) |
| `upstream/` | `00_reference.sh` … `07_multiqc.sh`: one shell script per step, FASTQ → counts |
| `downstream/` | `00_upstream_summary.R` … `12_ReactomePA.R`: one R script per results folder |
| `slides/` | the presentation source (`rnaseq_intro.qmd`, Quarto/reveal.js) and `render.sh` |
| `results/` | the outputs (figures and tables), one folder per step; see [What you get](#what-you-get) |
| `docs/` | the rendered slides |
| `data/` | reads, reference genome, alignments; created by the pipeline, **not** in git (~11 GB) |

## Quick start

```bash
git clone https://github.com/mortezasabri/rnaseq-workshop.git && cd rnaseq-workshop
micromamba env create -f environment.yml      # once, ~5 min (or: conda env create -f environment.yml)
micromamba activate rnaseq-workshop
./run_all.sh --dry-run                        # checks tools, R packages, settings, sample sheet
./run_all.sh                                  # everything: ~15 min, plus a one-time 6 GB download
```

Laptop only, or no Linux/Mac? Skip the upstream part and start from the published count table
(the same 8 samples at full depth, from the `airway` R package). Needs only R:

```bash
Rscript install.R                             # once, if you do not use the conda environment
COUNTS_SOURCE=airway ./run_all.sh downstream
COUNTS_SOURCE=airway ./run_all.sh slides
```

Each step can also run on its own, from the repository folder: `bash upstream/03_trim.sh`,
`Rscript downstream/04_go.R`, or line by line in RStudio (set the working directory to the
repository folder first). Each R step reads what the step before it saved in
`results/downstream/objects/`.

## Requirements

| | |
|---|---|
| System | Linux or macOS for the upstream part (Apple Silicon: if a tool will not install, create the environment with `CONDA_SUBDIR=osx-64`). The downstream part runs anywhere R runs. |
| Software | conda, mamba or micromamba; everything else comes from `environment.yml` (R 4.5, Bioconductor 3.22, HISAT2, Salmon, …) |
| Disk | ~11 GB: human genome index 4.2 GB, Salmon index 1.6 GB, reads + alignments ~5 GB |
| Memory | 8 GB (HISAT2 holds the human genome index in memory) |
| Internet | Downloads ~6 GB on the first run (genome index, annotation, reads). KEGG, WikiPathways and Europe PMC are queried online on every run. |
| Time | Measured on 8 threads: upstream ~10 min (4 of them downloading reads), downstream ~5 min, slides seconds. The first run also downloads the 6 GB reference (~15 min at 6 MB/s). |

## What you get

`results/upstream/`: `multiqc/multiqc_report.html` (all QC in one report), `featureCounts/gene_counts.tsv`
(the count matrix), `salmon/gene_counts_salmon.tsv`, and the per-tool logs.

`results/downstream/`: one folder per step. Every folder holds PNG figures and the matching tables (TSV).

The repository already contains the count tables and every downstream figure and table from one run,
so you can browse the results before running anything. The QC report and the logs appear once you run
the upstream part.

| Folder | Content |
|---|---|
| `0_upstream` | read funnel (trim → align → count), read coverage of one gene, featureCounts vs. Salmon, pipeline vs. published counts |
| `1_eda` | library sizes, normalization, PCA, sample-to-sample distances |
| `2_dea` | DESeq2 results (all genes, significant genes, summary), DESeq2 vs. edgeR, subset vs. full-depth data |
| `3_dea_vis` | MA plot, volcano plot, heatmap of the top 30 genes, top genes per donor, p-value histogram, plain vs. shrunken fold changes |
| `4_go` | GO over-representation: dot plots (all ontologies, BP), bar plot, gene-concept network |
| `5_go_ora_clust` | up vs. down compared, redundant terms removed, enrichment map, term tree |
| `6_go_gse_clust` | GO gene set enrichment analysis (GSEA): dot plot, running-score plot, ridge plot, map, tree |
| `7_enriched_path` | KEGG over-representation and GSEA, plus the top pathway as a colored KEGG map |
| `8_do_ora` | Disease Ontology over-representation |
| `9_do_gse` | Disease Ontology GSEA |
| `10_pubmed` | publication trends (Europe PMC) for the top disease terms and top genes |
| `11_wikiPath` | WikiPathways over-representation and GSEA |
| `12_ReactomePA` | Reactome over-representation and GSEA, one pathway as a gene network |

A step that finds no significant term writes `<table>_NONE.txt` instead of an empty plot.

## Input: the sample sheet

`config/samples.tsv` lists the 8 airway samples. For your own data, copy the template and fill it in:

```bash
cp config/samples.template.tsv config/my_samples.tsv     # then edit it
SAMPLES=config/my_samples.tsv ./run_all.sh --dry-run
SAMPLES=config/my_samples.tsv ./run_all.sh
```

| Column | Meaning |
|---|---|
| `sample` | unique short name (letters, digits, `_`, `-`) |
| `condition` | the group; exactly two values, one of them `REFERENCE_LEVEL` (default `untreated`) |
| `batch` | donor, cell line, or sequencing batch; added to the model as `~ batch + condition`. `-` if none |
| `fastq_1` | read 1: local path (`.fastq` or `.fastq.gz`) or http(s)/ftp URL |
| `fastq_2` | read 2 for paired-end data, `-` for single-end |

The reference is human (GRCh38, Ensembl). Other species need other URLs in `config/config.sh`
and a different annotation package in the R scripts (`org.Hs.eg.db`). Set `N_READS=0` to use every read.

## Options

All settings are in `config/config.sh`. Override any of them for one run on the command line, for
example `THREADS=16 N_READS=500000 ./run_all.sh`.

| Setting | Default | Meaning |
|---|---|---|
| `SAMPLES` | `config/samples.tsv` | sample sheet |
| `N_READS` | `2000000` | reads (pairs) taken from the start of each FASTQ; `0` = all |
| `THREADS` | `8` | CPU threads |
| `STRANDED` | `0` | library strandedness for featureCounts (0 / 1 / 2); Salmon reports the detected type |
| `TRIM_WINDOW`, `TRIM_MINLEN` | `4:20`, `25` | Trimmomatic quality window and minimum read length |
| `RUN_SALMON` | `yes` | also quantify with Salmon (alignment-free) |
| `COUNTS_SOURCE` | `pipeline` | `pipeline` = our featureCounts table, `airway` = published full-depth table |
| `REFERENCE_LEVEL` | `untreated` | the baseline condition |
| `PADJ_CUTOFF`, `LFC_CUTOFF` | `0.05`, `1` | significance: adjusted p-value and \|log2 fold change\| |
| `SHOWCASE_GENE` | `CRISPLD2` | gene shown as read coverage (the gene highlighted in the airway paper) |
| `SEED` | `42` | random seed for GSEA, so reruns give identical results |

## What success and failure look like

`run_all.sh` checks everything first (tools, R packages, sample sheet, input files) and runs nothing
if a check fails. Then it prints one line per step:

```
[21:31:40] downstream  04_go                    OK      (0.6 min)
```

| Exit code | Meaning | What to do |
|---|---|---|
| 0 | all steps finished | results are in `results/`, slides in `docs/` |
| 1 | setup problem; nothing was run | follow the listed fixes; usually `micromamba activate rnaseq-workshop` |
| 2 | a step failed | the last lines of its log are printed; full log in `logs/<part>_<step>.log` |

- The upstream part stops at the first failed step, since each step needs the previous one.
- In the downstream part a failed enrichment step (for example KEGG offline) does not stop the
  others; the run ends with exit code 2 and names the failed steps.
- One line per step is also appended to `logs/run_all.log`.
- Upstream steps skip work that is already done: rerunning after a failure continues where it stopped.
  Each downstream step empties and rewrites its own results folder.

## Data

All data is public (ENA run accessions SRR1039508–SRR1039521) and contains nothing
patient-identifiable. Local FASTQ files listed in a sample sheet are read in place and never uploaded
anywhere; the only network traffic is downloading the reference, the airway reads, and the online
enrichment databases (gene IDs only). `data/` and the intermediate R objects are not tracked by git.

## Suggested 25-minute talk

| Time | Slides | Content |
|---|---|---|
| 3 min | 1–3 | The question, the experiment, the road map |
| 8 min | 4–10 | Upstream: FASTQ and quality scores, download and QC, trim/align/count, coverage, count matrix, Salmon, which tool when |
| 8 min | 11–17 | Normalization, PCA, the DESeq2 model, MA/volcano, gene by gene, is a 2 M-read subset enough? |
| 6 min | 18–25 | ORA vs. GSEA, GO, term clusters, GSEA, KEGG map, Reactome/WikiPathways, diseases and literature |
| 1 min | 26–27 | Take-home messages, how to run it |

More figures are in the backup slides at the end.

## References

- Himes BE *et al.* RNA-Seq transcriptome profiling identifies CRISPLD2 as a glucocorticoid responsive gene
  that modulates cytokine function in airway smooth muscle cells. *PLoS One* 2014;9:e99625.
- Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with
  DESeq2. *Genome Biol* 2014;15:550. Workflow: <https://bioconductor.org/packages/rnaseqGene>
- Yu G. Thirteen years of clusterProfiler. *The Innovation* 2024;5:100722.
  Book: <https://yulab-smu.top/biomedical-knowledge-mining-book/>
- Kim D *et al.* (HISAT2) *Nat Biotechnol* 2019;37:907–915. Liao Y *et al.* (featureCounts)
  *Bioinformatics* 2014;30:923–930. Patro R *et al.* (Salmon) *Nat Methods* 2017;14:417–419.
  Dobin A *et al.* (STAR) *Bioinformatics* 2013;29:15–21. Bray NL *et al.* (kallisto) *Nat Biotechnol*
  2016;34:525–527.
- Kanehisa M, Goto S. KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Res* 2000;28:27–30.
  Kanehisa M *et al.* KEGG: biological systems database as a model of the real world. *Nucleic Acids Res*
  2025;53:D672–D677.
- Milacic M *et al.* The Reactome Pathway Knowledgebase 2024. *Nucleic Acids Res* 2024;52:D672–D678.
  Yu G, He QY. ReactomePA: an R/Bioconductor package for reactome pathway analysis and visualization.
  *Mol Biosyst* 2016;12:477–479.
- Agrawal A *et al.* WikiPathways 2024: next generation pathway database. *Nucleic Acids Res*
  2024;52:D679–D689.
- Yu G, Wang LG, Yan GR, He QY. DOSE: an R/Bioconductor package for disease ontology semantic and
  enrichment analysis. *Bioinformatics* 2015;31:608–609.

## License

- Code (the pipeline, shell and R scripts): MIT, see [`LICENSE`](LICENSE).
- Slides (`slides/`, `docs/`) and the figures and tables in `results/`: Creative Commons Attribution 4.0
  ([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)). Reuse and adapt them freely, with credit.
- The airway data belongs to its authors (Himes *et al.* 2014, GEO GSE52778).
