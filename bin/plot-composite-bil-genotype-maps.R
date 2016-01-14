# Import functions for clustering and plotting samples based on bin genotypes
# From: https://github.com/mfcovington/detect-boundaries (commit b5de0b6)
source("plot-composite-map.R")

par1 <- "M82"
par2 <- "PEN"

figure.num <- 4
figure.dir <- file.path("figures", figure.num)
if (!dir.exists(figure.dir))
  dir.create(figure.dir)


# Will plot bin genotypes using either physical or genetic distance
bins.physical.file <- "data/bins/bin-genotypes.BILs.2014-12-07.imputed-NAs.merged-like"
bins.genetic.file <- "data/bins/bin-genotypes.BILs.2014-12-07.imputed-NAs.merged-like.genetic-distance"

bins.physical <- read.table(bins.physical.file, header = T, sep = "\t")
bins.genetic <- read.table(bins.genetic.file, header = T, sep = "\t")


# Use abbreviated chromosome IDs (e.g. 'ch01' instead of 'SL2.40ch01')
bins.physical$chr <- sub(".*(ch\\d+)", "\\1", bins.physical$chr)
bins.genetic$chr <- sub(".*(ch\\d+)", "\\1", bins.genetic$chr)


# Cluster based on genotypes across entire genome
order <- GetOrderOfSamplesClusteredByGenotype(bins.physical,
                                              par1 = par1, par2 = par2)

bins.physical.m <- ClusterAndMeltBinGenotypes(bins.physical, order,
                                              par1 = par1, par2 = par2)
bins.genetic.m <- ClusterAndMeltBinGenotypes(bins.genetic, order,
                                             par1 = par1, par2 = par2)

plot.file <- file.path(figure.dir, "composite-map.physical.png")
PlotCompositeMap(bins.physical.m, par1 = par1, par2 = par2, col1 = "magenta",
                 col2 = "green", plot.file = plot.file,
                 save = TRUE, plot=FALSE, chr.text.size = 12, width = 10,
                 height = 7.5)

plot.file <- file.path(figure.dir, "composite-map.genetic.png")
PlotCompositeMap(bins.genetic.m, par1 = par1, par2 = par2, col1 = "magenta",
                 col2 = "green", plot.file = plot.file,
                 save = TRUE, plot=FALSE, chr.text.size = 12, width = 10,
                 height = 7.5, genetic.distance = TRUE)


# Cluster by chromosome after removing samples without an introgression
bins.physical.m <- data.frame()
bins.genetic.m <- data.frame()

for (chromosome in unique(bins.physical$chr)) {
  bins.physical.chr <- SubsetByChrAndHasIntrogression(bins.physical, chromosome)
  bins.genetic.chr <- SubsetByChrAndHasIntrogression(bins.genetic, chromosome)

  order <- GetOrderOfSamplesClusteredByGenotype(bins.physical.chr,
                                                par1 = par1, par2 = par2)

  bins.physical.chr.m <- ClusterAndMeltBinGenotypes(bins.physical.chr, order,
                                                    par1 = par1, par2 = par2)
  bins.genetic.chr.m <- ClusterAndMeltBinGenotypes(bins.genetic.chr, order,
                                                   par1 = par1, par2 = par2)

  bins.physical.m <- rbind(bins.physical.m, bins.physical.chr.m)
  bins.genetic.m <- rbind(bins.genetic.m, bins.genetic.chr.m)
}

plot.file <- file.path(figure.dir, "composite-map.physical.cluster-by-chr.png")
PlotCompositeMap(bins.physical.m, stacked.chromosomes = TRUE,
                 par1 = par1, par2 = par2, col1 = "magenta", col2 = "green",
                 plot.file = plot.file,
                 save = TRUE, plot = FALSE, chr.text.size = 12, width = 7.5,
                 height = 10)

plot.file <- file.path(figure.dir, "composite-map.genetic.cluster-by-chr.png")
PlotCompositeMap(bins.genetic.m, stacked.chromosomes = TRUE,
                 par1 = par1, par2 = par2, col1 = "magenta", col2 = "green",
                 plot.file = plot.file,
                 save = TRUE, plot = FALSE, chr.text.size = 12, width = 7.5,
                 height = 10, genetic.distance = TRUE)

