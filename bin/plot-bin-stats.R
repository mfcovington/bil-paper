# Import functions for calculating and plotting bin stats
# From: https://github.com/mfcovington/detect-boundaries (commit b5de0b6)
source("plot-bin-stats.R")

boundaries.dir <- "data/boundaries/"

par1 <- "M82"
par2 <- "PEN"

figure.num <- 3
figure.dir <- file.path("figures", figure.num)
if (!dir.exists(figure.dir))
  dir.create(figure.dir)


# Accumulate counts and lengths data for introgressions
CountAndMeasureIntrogressions(boundaries.dir, par1 = par1, par2 = par2)


# Plot number of introgressions per sample
plot.file <- file.path(figure.dir, "introgressions-per-sample.png")
PlotIntrogressionsPerSample(counts.df, save = TRUE, plot = FALSE,
                            plot.file = plot.file,
                            ggtitle = "Introgression Frequency per BIL",
                            ylab = "Number of BILs", width = 5, height = 7.5)


# Plot percent Solanum pennelli in genome
plot.file <- file.path(figure.dir, "percent-pennelli.png")
PlotPercentIntrogressed(lengths.df, save = TRUE, plot = FALSE,
                        plot.file = plot.file,
                        ggtitle = "Introgression Proportion per BIL",
                        xlab = "Percent of S. pennellii in Genome", ylab = "Number of BILs",
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

plot.file <- file.path(figure.dir, "distribution-of-introgressions.physical.png")
PlotDistributionOfIntrogressions(
      bins.physical, par1 = par1, par2 = par2, color.introgression = "green",
      plot.file = plot.file,
      ylab = "Number of BILs with Introgression",
      save = TRUE, plot = FALSE, width = 7.5, height = 10)

plot.file <- file.path(figure.dir, "distribution-of-introgressions.genetic.png")
PlotDistributionOfIntrogressions(
      bins.genetic, par1 = par1, par2 = par2, color.introgression = "green",
      plot.file = plot.file,
      ylab = "Number of BILs with Introgression",
      save = TRUE, plot = FALSE, width = 7.5, height = 10,
      genetic.distance = TRUE)


# Plot bins per chromosome
plot.file <- file.path(figure.dir, "bins-per-chromosome.png")
PlotBinsPerChromosome(
      bins.physical, plot.file = plot.file,
      save = TRUE, plot = FALSE, width = 7.5, height = 10)
