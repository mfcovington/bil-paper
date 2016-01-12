library(plyr)
library(ggplot2)

boundaries.dir <- "data/boundaries/"

par1_id <- "M82"
par2_id <- "PEN"


# Accumulate counts and lengths data for introgressions
counts.df <- data.frame(id = character(),
                        par1 = integer(),
                        het = integer(),
                        par2 = integer())
lengths.df <- counts.df

filelist <- list.files(boundaries.dir, pattern = ".+\\.boundaries")

for (filename in filelist) {
  id <- sub("(^.+)\\.boundaries$", "\\1", filename)

  df <- read.table(paste0(boundaries.dir, filename),
                   colClasses = c("factor", rep("numeric", 2), "factor"))
  colnames(df) <- c("chr", "start", "end", "genotype")
  df$length <- df$end - df$start + 1

  counts <- count(df, "genotype")
  par1_counts <- max(counts$freq[counts$genotype == par1_id], 0)
  het_counts  <- max(counts$freq[counts$genotype == "HET"], 0)
  par2_counts <- max(counts$freq[counts$genotype == par2_id], 0)

  counts.df <- rbind(counts.df, data.frame(id = id,
                                           par1 = par1_counts,
                                           het = het_counts,
                                           par2 = par2_counts))

  lengths <- aggregate(length ~ genotype, data = df, sum)
  par1_lengths <- max(lengths$length[lengths$genotype == par1_id], 0)
  het_lengths  <- max(lengths$length[lengths$genotype == "HET"], 0)
  par2_lengths <- max(lengths$length[lengths$genotype == par2_id], 0)

  lengths.df <- rbind(lengths.df, data.frame(id = id,
                                             par1 = par1_lengths,
                                             het = het_lengths,
                                             par2 = par2_lengths))
}


# Plot number of introgressions per sample
max.count.introgression.combined <- max(counts.df$het + counts.df$par2)

ggplot(data = counts.df) +
  geom_histogram(aes(x = het + par2),
                 binwidth = 1,
                 origin = 0.5,
                 color = 'black',
                 fill = 'skyblue') +
  xlim(0, max.count.introgression.combined + 0.5) +
  xlab("# of introgressions per sample") +
  ylab("# of BILs")

ggsave("figures/introgressions-per-sample.png", width = 5, height = 7.5)


ggplot(data=lengths.df) +
  geom_density(aes(x = par1)) +
  geom_density(aes(x = par2))

ggplot(data=lengths.df) +
  geom_histogram(aes(x = 100 * (par2 + het) / par1), binwidth = 2) +
  xlab("% of genome comprised of introgressions") +
  ylab("# of BILs")

# BAD DENOMINATOR
ggplot(data=lengths.df) +
  geom_histogram(aes(x = 100 * (par2 + het) / par1), binwidth = 0.15, color = 'black', fill = 'skyblue') +
  xlab("% of S. penn. in genome") +
  ylab("# of BILs") +
  scale_x_log10()

# THIS ONE?
ggplot(data=lengths.df) +
  geom_histogram(aes(x = 100 * (par2 + het) / (par1 + het + par2)), binwidth = 0.15, color = 'black', fill = 'skyblue') +
  xlab("% of S. penn. in genome") +
  ylab("# of BILs") +
  scale_x_log10()

# WITHOUT HET
ggplot(data=lengths.df) +
  geom_histogram(aes(x = 100 * par2 / (par1 + het + par2)), binwidth = 0.15, color = 'black', fill = 'skyblue') +
  xlab("% of S. penn. in genome") +
  ylab("# of BILs") +
  scale_x_log10()


#####

bin.geno.file <- "data/bins/bin-genotypes.BILs.2014-12-07.imputed-NAs.merged-like"
bin.geno.df <- read.table(bin.geno.file, header = TRUE, sep = "\t")

bin.geno.df$par1 <- apply(bin.geno.df, 1, function(line) sum(line == par1_id))
bin.geno.df$het <- apply(bin.geno.df, 1, function(line) sum(line == "HET"))
bin.geno.df$par2 <- apply(bin.geno.df, 1, function(line) sum(line == par2_id))

ggplot(bin.geno.df, aes(xmin = bin.start, xmax = bin.end),) +
  geom_rect(aes(ymin = 0, ymax = par2), fill = 'green') +
  geom_rect(aes(ymin = par2, ymax = het + par2), fill = 'black') +
  facet_grid(chr ~ .)


# added due to broken feature in ggplot 0.9.1:
# is it still broken/required?
max_pos <- max(bin.geno.df$bin.end)
temp_max <- floor(max_pos / 10000000)
max_pos_fixed <- temp_max * 10000000
label_max <- temp_max * 10

offset <- 1.05
max_count <-  max(bin.geno.df$het + bin.geno.df$par2)
max_lab <- (max_count %/% 10) * 10

# THIS ONE?
ggplot(bin.geno.df, aes(xmin = bin.start, xmax = bin.end),) +
  # geom_rect(aes(ymin = par2 + het, ymax = max_count), fill = 'magenta') +
  # geom_rect(aes(ymin = par2 + 5, ymax = het + par2 + 5), fill = 'black') +
  # geom_rect(aes(ymin = 0, ymax = het + par2), fill = 'black') +
  geom_rect(aes(ymin = 0, ymax = par2), fill = 'green') +
  geom_rect(aes(ymin = par2, ymax = het + par2), fill = 'black') +
  facet_grid(chr ~ .) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  ggtitle('Distribution of Introgressions Across Bins') +
  scale_x_continuous(
    'Bin position on chromosome (Mb)',
    breaks = seq(0, max_pos_fixed, 10000000),
    labels = seq(0, label_max,     10)
  ) +
  scale_y_continuous(
    '# of BILs with introgression',
    breaks = c(0, max_lab / 2, max_lab),
    labels = c(0, max_lab / 2, max_lab),
    limits = c(0, offset * max_count)
  ) +
  theme(strip.text.y = element_text(size = 5))
  # theme(strip.text.y = element_text(size = 5, angle = 90))
