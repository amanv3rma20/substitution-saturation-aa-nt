# substitution-saturation-aa-nt
R scripts for assessing nucleotide and amino acid substitution saturation. Compares uncorrected p-distances with model-corrected distances across codon positions and phylogenetic partitions, with length-weighted aggregate analyses and regression-based saturation plots.
# Substitution Saturation Analysis: Nucleotide and Amino Acid Datasets

R scripts for evaluating substitution saturation in nucleotide and amino acid sequence datasets used for phylogenetic analyses.

The scripts compare **uncorrected p-distances** with **model-corrected evolutionary distances** using linear regression to assess the extent to which multiple substitutions may obscure phylogenetic signal. Both nucleotide and amino acid datasets can be analyzed at the partition level, with additional length-weighted analyses for assessing saturation across combined datasets.

## Scripts

### `nt_saturation_by_partition.R`

Performs nucleotide substitution saturation analysis across phylogenetic partitions and codon positions.

The script:

* Calculates uncorrected nucleotide p-distances.
* Calculates TN93 model-corrected nucleotide distances.
* Separately analyzes:

  * All nucleotide sites
  * First codon positions
  * Second codon positions
  * Third codon positions
* Generates saturation plots and regression statistics for individual partitions.
* Performs length-weighted analyses across the combined dataset.
* Performs low-p saturation analyses and meta-analysis of partition-specific slopes.
* Produces PNG and SVG figures.

### `aa_saturation_partitioned.R`

Performs amino acid substitution saturation analysis for partitioned protein alignments.

The script:

* Calculates uncorrected amino acid p-distances.
* Calculates model-corrected amino acid distances.
* Supports multiple amino acid substitution models.
* Allows substitution models to be specified globally or independently for each partition.
* Generates partition-specific saturation plots.
* Calculates regression slopes and R² values for each partition.
* Produces pooled saturation statistics across all partitions.
* Produces PNG and SVG figures.

### `summarize_partition_slopes.R`

Aggregates the amino acid partition-level results produced by `aa_saturation_partitioned.R`.

The script:

* Reads partition-specific pairwise distance tables.
* Applies optional R² and partition-length filters.
* Calculates length-weighted mean distances for each taxon pair.
* Performs an overall regression between uncorrected and model-corrected distances.
* Generates a length-weighted aggregate saturation plot.

## Requirements

The scripts require R and the following packages:

* `ape`
* `phangorn`
* `tools`
* `ggplot2` *(optional; used for summary plots when available)*

Install the required packages in R with:

```r
install.packages(c("ape", "phangorn", "ggplot2"))
```

## Usage

### Nucleotide saturation

```bash
Rscript nt_saturation_by_partition.R alignment.fasta output_directory partition.txt
```

Arguments:

1. `alignment.fasta` — nucleotide alignment in FASTA format.
2. `output_directory` — directory for analysis outputs.
3. `partition.txt` — partition file defining the genomic regions or marker genes.

### Amino acid saturation

```bash
Rscript aa_saturation_partitioned.R alignment.fasta output_directory partition.txt [MODEL|AUTO]
```

The optional `MODEL` argument specifies a single amino acid substitution model for all partitions.

Alternatively, use:

```bash
Rscript aa_saturation_partitioned.R alignment.fasta output_directory partition.txt AUTO
```

to use the model specified for each partition in the partition file.

### Summarize amino acid partitions

```bash
Rscript summarize_partition_slopes.R pairs_directory partition_slopes.tsv output_directory [r2_min] [len_min]
```

Optional filters can be used to retain only partitions meeting a minimum R² value and/or minimum partition length.

## Interpretation

For both nucleotide and amino acid analyses, the **x-axis represents the uncorrected p-distance**, while the **y-axis represents the model-corrected evolutionary distance**.

Uncorrected p-distance represents the observed proportion of different sites between two sequences. Model-corrected distance accounts for substitutions that may have occurred multiple times at the same site.

A linear relationship with a high R² indicates that observed sequence divergence remains a good predictor of model-corrected evolutionary divergence. Deviations from linearity, particularly flattening at higher p-distances, can indicate increasing substitution saturation.

For nucleotide datasets, comparisons among codon positions can be used to identify position-specific saturation, particularly at third codon positions. Amino acid analyses provide a complementary assessment of saturation at the protein level.

## Output

The scripts generate:

* Pairwise distance tables (`.tsv`)
* Partition-specific regression statistics
* Aggregate regression statistics
* Saturation plots (`.png`)
* Vector saturation plots (`.svg`)
* Partition-level summary tables
* Low-p meta-analysis results

SVG files provide scalable vector graphics suitable for editing and preparation of publication-quality figures.

## Citation

If you use these scripts in your research, please cite the associated publication or repository:

**Verma, A.** *[Publication details to be added]*

Repository:
https://github.com/amanv3rma20/substitution-saturation-aa-nt
