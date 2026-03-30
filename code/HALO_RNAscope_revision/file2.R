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
  dt[, x_center := (XMin + XMax) / 2]
  
  # relative depth within each image: 0 = top, 1 = bottom
  xr <- range(dt$x_center, na.rm = TRUE)
  dt[, x_rel := (x_center - xr[1]) / (xr[2] - xr[1])]
  
  dt
})

combined_dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)

# keep only needed columns and convert to long format
plot_dt <- melt(
  combined_dt,
  id.vars = c("brain", "region", "file", "Object Id", "Analysis Region", "x_center", "x_rel"),
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
plot_dt <- plot_dt[positive == 1]
nbins <- 1000
plot_dt[, x_bin := cut(x_rel, breaks = seq(0, 1, length.out = nbins + 1), include.lowest = TRUE)]

# get numeric bin center
bin_breaks <- seq(0, 1, length.out = nbins + 1)
bin_centers <- (bin_breaks[-1] + bin_breaks[-length(bin_breaks)]) / 2

bin_map <- data.table(
  x_bin = levels(plot_dt$x_bin),
  x_bin_center = bin_centers
)

# count cells per region x marker x bin
count_dt <- plot_dt[, .N, by = .(region, celltype, x_bin)]
count_dt <- merge(count_dt, bin_map, by = "x_bin", all.x = TRUE)

# normalize within each region and marker so curves are comparable
count_dt[, density_scaled := N / sum(N), by = .(region, celltype)]
count_dt[, frac := N / sum(N), by = .(region, x_bin)]
count_dt[, x_bin := factor(x_bin, levels = levels(plot_dt$x_bin))]

# colors requested
celltype_cols <- c(
  "ADCYAP1+" = "green",
  "RORB+"    = "blue",
  "PCP4+"    = "pink",
  "AQP4+"    = "gold"
)

p_all <- ggplot(count_dt,
                aes(x = x_bin, y = frac, fill = celltype)) +
  geom_col(width = 1) +
  scale_fill_manual(values = celltype_cols) +
  facet_wrap(~region, ncol = 2) +
  labs(
    x = "X position (binned)",
    y = "Fraction of cells",
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p_all

p_density <- ggplot(plot_dt, aes(x = x_rel, fill = celltype, color = celltype)) +
  # alpha makes the overlapping areas visible
  geom_density(alpha = 0.3, linewidth = 0.8) + 
  scale_fill_manual(values = celltype_cols) +
  scale_color_manual(values = celltype_cols) +
  facet_wrap(~region, ncol = 2) +
  labs(
    x = "Relative X Position",
    y = "Density",
    fill = "Cell Type",
    color = "Cell Type"
  ) +
  theme_classic(base_size = 14)

p_density

p_line <- ggplot(count_dt, 
                 aes(x = x_bin_center, y = density_scaled, color = celltype)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = celltype_cols) +
  facet_wrap(~region, ncol = 2) +
  labs(
    x = "X position (center of bin)",
    y = "Scaled Density (N / Total N per type)",
    color = NULL
  ) +
  theme_classic(base_size = 14)

p_line

# You can adjust this value and re-run the binning logic above
# nbins <- 100 

ggplot(count_dt, aes(x = x_bin_center, y = density_scaled, color = celltype, fill = celltype)) +
  # 1. Semi-transparent area fill (position = "identity" is key to allow overlap)
  #geom_area(position = "identity", alpha = 0.2, color = NA) +
  stat_smooth(geom = "area", method = "loess", span = 0.1, alpha = 0.2, position = "identity", color = NA) +
  
  # 2. Sharp line for the distribution outline
  #geom_line(linewidth = 0.8) +
  stat_smooth(geom = "line", method = "loess", span = 0.1, linewidth = 0.8)+
  # 3. Facet by brain region
  facet_wrap(~region, ncol = 1, scales = "free_y") +
  
  # 4. Styling and Scales
  scale_fill_manual(values = celltype_cols) +
  scale_color_manual(values = celltype_cols) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  
  # 5. Labels
  labs(
    x = "Relative Depth (0 = Top, 1 = Bottom)",
    y = "Density (Fraction of cell type)",
    title = "Spatial Density of Cell Types across dACC and dlPFC",
    subtitle = paste0("Binned into ", nbins, " segments along the X-axis"),
    fill = "Cell Type",
    color = "Cell Type"
  ) +
  
  # 6. Theme Refinement
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 16),
    legend.position = "top",
    panel.spacing = unit(1, "lines")
  )
  
  
# dACC plot
p_dacc <- ggplot(count_dt[region == "dACC"],
                 aes(x = x_bin_center, y = N, color = celltype)) +
  geom_line(linewidth = 1.2, alpha = 0.2) +
  scale_color_manual(values = celltype_cols) +
  scale_y_reverse() +
  labs(
    title = "dACC",
    y = "Relative density",
    x = "Cortical depth (top to bottom)",
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