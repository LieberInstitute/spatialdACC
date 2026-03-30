library(data.table)
library(ggplot2)

indir <- "/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/raw-data/HALO_RNAscope_revision"

files <- list.files(indir, pattern = "\\.csv$", full.names = TRUE)

# read all files
dt_list <- lapply(files, function(f) {
  dt <- fread(f)
  fname <- basename(f)
  
  brain <- sub("^Br([0-9]+)_.*", "Br\\1", fname)
  region <- ifelse(grepl("_dACC_", fname, ignore.case = TRUE), "dACC",
                   ifelse(grepl("_dlPFC_", fname, ignore.case = TRUE), "dlPFC", NA))
  
  dt[, file := fname]
  dt[, brain := brain]
  dt[, region := region]
  
  # cell y midpoint
  dt[, y_center := (YMin + YMax) / 2]
  
  # relative depth within each image: 0 = top, 1 = bottom
  yr <- range(dt$y_center, na.rm = TRUE)
  dt[, y_rel := (y_center - yr[1]) / (yr[2] - yr[1])]
  
  dt
})

combined_dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)

# keep only needed columns and convert to long format
plot_dt <- melt(
  combined_dt,
  id.vars = c("brain", "region", "file", "Object Id", "Analysis Region", "y_center", "y_rel"),
  measure.vars = c("ADCYAP1+", "RORB+", "PCP4+", "AQP4+"),
  variable.name = "celltype",
  value.name = "positive"
)

# clean marker names
#plot_dt[, celltype := factor(celltype,
#                             levels = c("ADCYAP1.", "RORB.", "PCP4.", "AQP4."),
#                             labels = c("ADCYAP1+", "RORB+", "PCP4+", "AQP4+"))]

# keep only positive cells
plot_dt_dACC <- plot_dt[positive == 1 & region =="dACC"]
plot_dt_dlPFC <- plot_dt[positive == 1 & region =="dlPFC"]

# make y bins
nbins <- 2732
plot_dt[, y_bin := cut(y_rel, breaks = seq(0, 1, length.out = nbins + 1), include.lowest = TRUE)]

# get numeric bin center
bin_breaks <- seq(0, 1, length.out = nbins + 1)
bin_centers <- (bin_breaks[-1] + bin_breaks[-length(bin_breaks)]) / 2

bin_map <- data.table(
  y_bin = levels(plot_dt$y_bin),
  y_bin_center = bin_centers
)

# count cells per region x marker x bin
count_dt <- plot_dt[, .N, by = .(region, celltype, y_bin)]
count_dt <- merge(count_dt, bin_map, by = "y_bin", all.x = TRUE)

# normalize within each region and marker so curves are comparable
count_dt[, density_scaled := N / sum(N), by = .(region, celltype)]

# colors requested
celltype_cols <- c(
  "ADCYAP1+" = "green",
  "RORB+"    = "blue",
  "PCP4+"    = "pink",
  "AQP4+"    = "gold"
)

# dACC plot
p_dacc <- ggplot(count_dt[region == "dACC"],
                 aes(x = y_bin_center, y = N, color = celltype)) +
  geom_line(linewidth = 1.2, alpha = 0.5) +
  scale_color_manual(values = celltype_cols) +
  scale_y_reverse() +
  labs(
    title = "dACC",
    x = "Relative density",
    y = "Cortical depth (top to bottom)",
    color = NULL
  ) +
  theme_classic(base_size = 14)

# dlPFC plot
p_dlpfc <- ggplot(count_dt[region == "dlPFC"],
                  aes(x = density_scaled, y = y_bin_center, color = celltype)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = celltype_cols) +
  scale_y_reverse() +
  labs(
    title = "dlPFC",
    x = "Relative density",
    y = "Cortical depth (top to bottom)",
    color = NULL
  ) +
  theme_classic(base_size = 14)

p_dacc
p_dlpfc