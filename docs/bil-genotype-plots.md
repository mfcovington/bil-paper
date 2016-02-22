<!-- MarkdownTOC -->

- [Plot BIL genotype and boundary data](#plot-bil-genotype-and-boundary-data)
    - [Plots from original round of sequencing](#plots-from-original-round-of-sequencing)

<!-- /MarkdownTOC -->


# Plot BIL genotype and boundary data

After identifying polymorphisms, genotyping individuals from the BIL population, and detecting introgression boundaries for each BIL, we generated plots that display both the genotype and the boundary data for each BIL.

Some of the samples from the original round of sequencing had very low expression or had some contamination from other samples. Where possible, these samples were resequenced.

## Plots from original round of sequencing

Plot BIL genotype and boundary data together using the R script `bin/plot-bil-genotypes-with-bins.R` as-is. This code requires that the working directory (`data`, in this case) contains the following input files (where `XXX` is the BIL sample ID:

- `boundaries/XXX.boundaries`
- `genotyped/XXX.genotyped.nr`

Plots are created in a sub-directory called `plots-with-boundaries`.
