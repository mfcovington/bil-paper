# Import functions for calculating and plotting bin stats
# From: https://github.com/mfcovington/detect-boundaries (commit d37d9a4)
source("plot-bin-stats.R")

boundaries.dir <- "data/boundaries/"

par1 <- "M82"
par2 <- "PEN"


# Accumulate counts and lengths data for introgressions
CountAndMeasureIntrogressions(boundaries.dir, par1 = par1, par2 = par2)


# Plot number of introgressions per sample
PlotIntrogressionsPerSample(counts.df, save = TRUE, plot = FALSE,
                            plot.file = "figures/introgressions-per-sample.png",
                            ylab = "# of BILs", width = 5, height = 7.5)


# Plot percent Solanum pennelli in genome
PlotPercentIntrogressed(lengths.df, save = TRUE, plot = FALSE,
                        plot.file = "figures/percent-pennelli.png",
                        xlab = "% of S. penn. in genome", ylab = "# of BILs",
                        width = 5, height = 7.5)


# Plot distribution of introgressions across bins
# using either physical or genetic distance
bins.physical.file <- "data/bins/bin-genotypes.BILs.2014-12-07.imputed-NAs.merged-like"
bins.genetic.file  <- "data/bins/bin-genotypes.BILs.2014-12-07.imputed-NAs.merged-like.genetic-distance"

bins.physical <- read.table(bins.physical.file, header = T, sep = "\t")
bins.genetic  <- read.table(bins.genetic.file, header = T, sep = "\t")

# Use abbreviated chromosome IDs (e.g. 'ch01' instead of 'SL2.40ch01')
bins.physical$chr <- sub(".*(ch\\d+)", "\\1", bins.physical$chr)
bins.genetic$chr <- sub(".*(ch\\d+)", "\\1", bins.genetic$chr)

PlotDistributionOfIntrogressions(
      bins.physical, par1 = par1, par2 = par2, color.introgression = "green",
      plot.file = "figures/distribution-of-introgressions.physical.png",
      save = TRUE, plot = FALSE, width = 7.5, height = 10)

PlotDistributionOfIntrogressions(
      bins.genetic, par1 = par1, par2 = par2, color.introgression = "green",
      plot.file = "figures/distribution-of-introgressions.genetic.png",
      save = TRUE, plot = FALSE, width = 7.5, height = 10,
      genetic.distance = TRUE)
