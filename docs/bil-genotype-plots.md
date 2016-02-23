<!-- MarkdownTOC -->

- [Plot BIL genotype and boundary data](#plot-bil-genotype-and-boundary-data)
    - [Plots from original round of sequencing](#plots-from-original-round-of-sequencing)
    - [Plots from resequenced BIL samples](#plots-from-resequenced-bil-samples)
    - [Aggregate plots from original and resequenced samples](#aggregate-plots-from-original-and-resequenced-samples)

<!-- /MarkdownTOC -->


# Plot BIL genotype and boundary data

After identifying polymorphisms, genotyping individuals from the BIL population, and detecting introgression boundaries for each BIL, we generated plots that display both the genotype and the boundary data for each BIL.

Some of the samples from the original round of sequencing had very low expression or had some contamination from other samples. Where possible, these samples were resequenced.

## Plots from original round of sequencing

Plot BIL genotype and boundary data together using the R script `bin/plot-bil-genotypes-with-bins.R` as-is. This code requires that the working directory (`data`, in this case) contains the following input files (where `XXX` is the BIL sample ID:

- `boundaries/XXX.boundaries`
- `genotyped/XXX.genotyped.nr`

Plots are created in a sub-directory called `plots-with-boundaries`.

## Plots from resequenced BIL samples

Plot BIL genotype and boundary data together using the R script `bin/plot-bil-genotypes-with-bins.R` with the following changes. The first two lines are replaced with:

```r
setwd("data/resequenced")
plot.title.supplement <- "resequenced"
```

Plots are created in a sub-directory called `resequenced/plots-with-boundaries`.


## Aggregate plots from original and resequenced samples

The final plots are located in a sub-directory called `plots-with-boundaries.FINAL`.

```sh
cd $OUT_DIR    # 'data', in this case

mkdir plots-with-boundaries.FINAL
cd plots-with-boundaries.FINAL/

# Make symbolic links to plots from the original round of sequencing
for PLOT in ../plots-with-boundaries/BIL_*pdf; do
    ln -s $PLOT
done

# Override symbolic links with plots from resequenced samples if resequencing resolves original issue
for ID in 015 020 029 032 035 038 043 044 047 095 097 110 120 133 148 155 177 182 183 185 187 242 243 291 294 304 310 349 350 359 389 396 405 406 408 429 432 473 480 497 510 540; do
    rm BIL_$ID.overlay.pdf
    ln -s ../resequenced/plots-with-boundaries/BIL_$ID.overlay.pdf
done

# Include plots from both original and resequenced samples where useful
# BIL_089, BIL_232, BIL_245: Both original and resequenced samples were contaminated with another BIL
# BIL_269: Original suggested large heterozygous region. Resequenced sample confirms this.
for ID in 089 232 245 269; do
    rm BIL_$ID.overlay.pdf
    ln -s ../plots-with-boundaries/BIL_$ID.overlay.pdf BIL_$ID.overlay.original.pdf
    ln -s ../resequenced/plots-with-boundaries/BIL_$ID.overlay.pdf BIL_$ID.overlay.resequenced.pdf
done
```

Plots with issues that were not resequenced, but for which there was enough information to infer introgression boundaries:

- BIL_328

    - No germination for resequencing
    - Contaminated with DNA from BIL_266

- BIL_476

    - No usable DNA for resequencing

- BIL_534

    - No germination for resequencing
    - Contaminated with DNA from BIL_102
