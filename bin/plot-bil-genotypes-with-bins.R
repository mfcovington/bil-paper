# Working directory must contain following:
# - boundaries/XXX.boundaries
# - genotyped/XXX.genotyped.nr
setwd("data")

plot.title.supplement <- ""    # Optional: Adds extra info inside parentheses

library(ggplot2)
library(reshape)

theme_set(theme_grey(base_size = 12))

par1 <- "M82"
par2 <- "PEN"

col_par1 <- "magenta"
col_par2 <- "green"
col_het  <- "black"
col_bound_line <- "white"
col_bound_fill <- "gray"

offset <- 1.05

boundary.files <- list.files(path = "boundaries", pattern = "\\.boundaries$")

output.dir <- "plots-with-boundaries"
if (!file.exists(output.dir)) dir.create(output.dir)

for (b.file in boundary.files) {
  id <- sub("(.+)\\.boundaries", "\\1", b.file)
  if (file.info(paste0("boundaries/", b.file))$size == 0) {
    print(paste("WARNING: NO BOUNDARIES FOR", id))
    next
  }

  geno.files <- list.files(path = "genotyped", pattern = paste0(id, "\\..+\\.genotyped\\.nr$"))
  if (length(geno.files) == 0) {
    print(paste("WARNING: NO GENOTYPED FILES FOR", id))
    next
  }

  print(paste("Plotting", id))

  bounds <- read.table(paste0("boundaries/", b.file))
  colnames(bounds) <- c("chr", "start", "end", "geno")

  geno_df <- data.frame()
  for (g.file in geno.files) {
      if (file.info(paste0("genotyped/", g.file))$size == 0) next
      geno_df <- rbind(geno_df,read.table(paste0("genotyped/", g.file)))
  }
  colnames(geno_df) <- c("chr","pos","par1", "par2", "cov")

  # Use shortened chromosome IDs
  bounds$chr  <- sub(".*(ch\\d+)", "\\1", bounds$chr)
  geno_df$chr <- sub(".*(ch\\d+)", "\\1", geno_df$chr)

  geno_df$cov.plot <- geno_df$cov
  geno_df$cov.plot[geno_df$par2 - geno_df$par1 < 0] <- -geno_df$cov.plot[geno_df$par2 - geno_df$par1 < 0]
  geno_df$cov.plot[geno_df$par1 + geno_df$par2 == 0] <- 0    # If all reads are non-parental, plot at y = 0

  # added due to broken feature in ggplot 0.9.1:
  # is it still broken/required?
  max_pos <- max(geno_df$pos)
  temp_max <- floor(max_pos / 10000000)
  max_pos_fixed <- temp_max * 10000000
  label_max <- temp_max * 10

  max <-  max(geno_df$cov.plot)
  min <- -min(geno_df$cov.plot)
  if(max < 0) max <- 0
  if(min < 0) min <- 0
  max_lab <- (max %/% 10) * 10
  min_lab <- (min %/% 10) * 10

  bounds$score[bounds$geno == par1]  <- -min
  bounds$score[bounds$geno == par2]  <-  max
  bounds$score[bounds$geno == "HET"] <-  0

  mbounds <- melt(bounds, id=c("chr", "score", "geno"))

  plot.title <- ifelse(nchar(plot.title.supplement) > 0,
                       paste0(id, " (", plot.title.supplement, ")"),
                       id)

  geno.plot <- ggplot() +
    geom_ribbon(data = mbounds, aes(x = value, ymin = 0, ymax = score),
                size = 0.4, color = col_bound_line, fill = col_bound_fill)
  if (sum(mbounds$geno == "HET") > 0) {
    geno.plot <- geno.plot +
      geom_rect(data = bounds[bounds$geno == "HET",], aes(xmin = start, xmax = end),
                ymin = -min, ymax = max, size = 0.4, color = col_bound_line, fill = col_bound_fill)
  }
  geno.plot <- geno.plot +
    geom_point(data = geno_df, shape = ".", aes(pos, cov.plot, color = (par2 - par1) / cov)) +
    facet_grid(chr ~ .) +
    scale_color_gradient2(
      low    = col_par1,
      mid    = col_het,
      high   = col_par2,
      limits = c(-1, 1),
      name   = substitute(frac(par2 - par1, total),
          list(par1 = par1, par2 = par2, total = "Total reads"))
    ) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      legend.title     = element_text(size = 10),
      strip.text.y     = element_text(size = 12)
    ) +
    ggtitle(plot.title) +
    scale_x_continuous(
      'Position on chromosome (Mb)',
      breaks = seq(0, max_pos_fixed, 10000000),
      labels = seq(0, label_max,     10)
    ) +
    scale_y_continuous(
      'Coverage',
      breaks = c(-min_lab, 0, max_lab),
      labels = c( min_lab, 0, max_lab),
      limits = c(-offset*min, offset*max)
    )

  ggsave(paste0(output.dir, "/", id, ".overlay.pdf"), plot = geno.plot, width = 10, height = 8)
}
