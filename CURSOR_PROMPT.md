# Cursor Implementation Prompt
# RNA-seq Differential Expression Pipeline — Project 1

---

## CONTEXT

You are helping implement a professional bioinformatics GitHub portfolio project.
This is a production-quality RNA-seq pipeline for a graduate student (M.S. Bioinformatics,
Johns Hopkins) targeting entry-level Bioinformatics Analyst / Computational Biologist roles.

The project owner has designed the architecture and provided starter code.
Your job is to **implement, complete, test, and polish** each file.

---

## REPOSITORY STRUCTURE (already scaffolded)

```
rna-seq-pipeline/
├── Snakefile
├── config/
│   ├── config.yaml
│   └── samples.tsv
├── rules/
│   ├── qc.smk
│   ├── trim.smk
│   ├── align.smk
│   ├── quantify.smk
│   └── deseq2.smk
├── scripts/
│   ├── deseq2_analysis.R
│   ├── pathway_enrichment.R
│   └── clean_count_matrix.py
├── envs/
│   ├── alignment.yaml
│   ├── quantification.yaml
│   └── r_analysis.yaml
├── notebooks/
│   └── exploratory_analysis.ipynb
└── README.md
```

---

## YOUR TASKS (in order)

### TASK 1 — Validate and fix the Snakefile

Open `Snakefile`. Verify:
- All `include:` paths resolve to existing `.smk` files
- The `rule all` input list correctly references output files from all rules
- `get_fastq_r1` and `get_fastq_r2` helper functions correctly index `SAMPLES`
- `SAMPLE_IDS` is derived correctly from `config/samples.tsv`

Fix any syntax errors. Add a `localrules: all` declaration.
Add a `wildcard_constraints` block to constrain `{sample}` to alphanumeric + underscore.

---

### TASK 2 — Complete rules/qc.smk

The rule `fastqc_raw` is complete. Verify:
- Output naming matches what MultiQC expects
- Log paths are consistent with other rules

Add a rule `samtools_stats` that runs `samtools stats` on each BAM and outputs
`{results}/qc/{sample}_stats.txt`. This feeds into MultiQC automatically.

---

### TASK 3 — Complete rules/align.smk

The STAR rules are scaffolded. Verify and fix:
- `star_genome_index` output sentinel is correct (Snakemake `directory()` rules require special handling)
- `star_align` uses `--outSAMtype BAM SortedByCoordinate` correctly (space in string)
- Add `--outBAMsortingBinsN 50` to prevent STAR temp file errors on large datasets
- Add `--limitBAMsortRAM 30000000000` (30 GB) for large genome alignment

---

### TASK 4 — Complete scripts/deseq2_analysis.R

The script is mostly complete. Fix these items:

1. The `lfcShrink` coef name (`"condition_tumor_vs_normal"`) must exactly match
   the DESeq2 results name. Add a line that prints `resultsNames(dds)` to the log
   so the user can verify the coef string is correct.

2. Add error handling: if `top50_genes` has fewer than 50 genes (small dataset),
   use however many are available without throwing an error.

3. Add a `sessionInfo()` call at the end of the log for reproducibility.

4. The `%||%` operator is not base R — replace it with an explicit `if (is.null(...))` check
   in the GSEA params, or add `"%||%" <- function(a, b) if (!is.null(a)) a else b` at the top.

---

### TASK 5 — Complete scripts/pathway_enrichment.R

Fix this syntax error:
```r
message(sprintf("Gene sets loaded: %d pathways", length(unique(msig_t2g$gs_name)))
```
There is a missing closing parenthesis. Fix it.

Also add:
- A check: if `n_sig == 0`, write a warning to the log and save an empty placeholder
  plot (a ggplot with `annotate("text", ...)` saying "No significant pathways found")
  so the pipeline does not fail.

---

### TASK 6 — Create notebooks/exploratory_analysis.ipynb

Create a Jupyter notebook with these sections:

**Section 1: Setup**
```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
```

**Section 2: Load DESeq2 results**
- Load `results/tables/deseq2_results.tsv`
- Print shape, head, value_counts of `direction` column

**Section 3: Load normalized counts**
- Load `results/tables/normalized_counts.tsv`
- Print shape

**Section 4: DEG summary stats**
- Total DEGs, upregulated count, downregulated count
- Top 10 most significant genes (by padj)
- Top 10 largest |log2FoldChange| genes

**Section 5: Volcano plot (matplotlib)**
- Recreate the volcano plot in Python/matplotlib as a demonstration
  that the candidate can work in both R and Python

**Section 6: Sample correlation heatmap**
- Compute pairwise Pearson correlation of normalized counts
- Plot as seaborn heatmap with sample labels and condition color bar

**Section 7: Gene of interest deep-dive**
- Pick 4 known breast cancer genes: `MKI67`, `ESR1`, `ERBB2`, `TP53`
- Plot their normalized expression as boxplots split by condition

Each section should have a Markdown cell explaining what is being done and why.

---

### TASK 7 — Add .gitignore

Create `.gitignore` with:
```
# Data files (large — do not commit to GitHub)
data/raw/
data/genome/
*.fastq.gz
*.bam
*.bai
*.fa
*.fa.gz
*.gtf
*.gtf.gz

# Snakemake internals
.snakemake/
__pycache__/
*.pyc

# R session artifacts
.Rhistory
.RData
*.Rproj

# Results (optional — uncomment if you want to track results)
# results/

# Conda environments (built locally)
envs/.conda/

# Jupyter checkpoints
.ipynb_checkpoints/
```

---

### TASK 8 — Validate pipeline syntax

Run this command from the project root and fix any errors until it exits cleanly:

```bash
snakemake --dry-run --cores 1 --use-conda 2>&1 | head -50
```

If you cannot run Snakemake in this environment, at minimum:
- Run `python -c "import ast; ast.parse(open('scripts/clean_count_matrix.py').read())"` to validate Python
- Run `Rscript --vanilla -e "parse(file='scripts/deseq2_analysis.R')"` to validate R syntax

---

### TASK 9 — Polish README.md

In the README:
1. Replace all `[Your Name]` and `YOUR_USERNAME` placeholders with a prompt:
   `<!-- TODO: Replace with your name and GitHub username -->`
2. Verify that the repository structure diagram matches the actual files on disk
3. Add a **Troubleshooting** section at the bottom with these common issues:
   - STAR OOM error (increase `--limitBAMsortRAM`)
   - featureCounts strand detection (how to determine strand from library prep protocol)
   - DESeq2 convergence warnings (what they mean, when to ignore them)
   - conda env solver slow (recommend `mamba` as drop-in replacement)

---

### TASK 10 — Final checklist before commit

Verify each item:
- [ ] `snakemake --dry-run` exits with 0 errors
- [ ] All `.smk` files have consistent `{results}` wildcard usage
- [ ] All log paths use `{results}/logs/{rule_name}/{sample}.log` pattern
- [ ] `config/samples.tsv` has correct tab-separated columns (not spaces)
- [ ] `.gitignore` exists and covers FASTQ, BAM, and genome files
- [ ] `README.md` has the Quick Start section and dataset instructions
- [ ] `notebooks/exploratory_analysis.ipynb` runs without errors on the results

---

## STYLE CONVENTIONS

- Python: PEP 8. Docstrings on all functions. Type hints preferred.
- R: tidyverse style. Comments on non-obvious lines. No `attach()`.
- Snakemake: snake_case rule names. All rules have `log:` and `conda:` directives.
- Commits: conventional commit format (`feat:`, `fix:`, `docs:`, `chore:`)

## DO NOT

- Do not download any real data — use placeholder paths
- Do not change the config structure (other scripts depend on it)
- Do not rename output files — README and rule targets depend on exact names
- Do not add heavy dependencies not in the envs/ YAML files

---

## WHEN YOU ARE DONE

The repository should pass `snakemake --dry-run` and contain no placeholder `TODO`
comments except in README.md where noted. All scripts should be syntactically valid.

---

## TASK 8 — Validation log (2025-06-01)

| Check | Result |
|---|---|
| `scripts/*.py` (`ast.parse`) | **PASS** — `clean_count_matrix.py` OK |
| `scripts/*.R` (`Rscript` parse) | **SKIPPED** — `Rscript` not on PATH in this environment |
| `snakemake --dry-run --cores 1 --use-conda` | **SKIPPED** — `snakemake` / `conda` not installed (pip install of snakemake 7.32 failed on Windows: `datrie` wheel build) |

**Manual review:** Snakefile `rule all` targets aligned with `rules/qc.smk` (R1/R2 FastQC HTML), `rules/deseq2.smk`, and `rules/quantify.smk`. MultiQC inputs fixed to use raw R1/R2 FastQC zips only (removed broken trimmed-zip dependency). R scripts edited per Tasks 2–3; parenthesis and guard blocks verified by inspection.

**Recommended on Linux/macOS or WSL:** `bash setup.sh && conda activate snakemake && snakemake --dry-run --cores 1 --use-conda`
