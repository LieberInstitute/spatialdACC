library(data.table)
library(ggplot2)

indir <- "/Users/madhavi.tippani/Desktop/HALO_dACC/ObjectData"

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
# Create spatial bins (1-100) for bar plots

# 1. Define the samples that need flipping (WM -> L1) 
# based on your table to bring them to L1 -> WM
combined_dt[, needs_flip := FALSE]

combined_dt[(brain == "Br6471" & region == "dACC") |
            (brain == "Br8492" & region == "dlPFC") |
            (brain == "Br2720" & region %in% c("dlPFC", "dACC")) |
            (brain == "Br6432" & region == "dACC") |
            (brain == "Br8667" & region == "dACC"), 
            needs_flip := TRUE]

# 2. Apply the flip: if it's WM -> L1, make it L1 -> WM
combined_dt[needs_flip == TRUE, x_rel := 1 - x_rel]

# keep only needed columns and convert to long format
plot_dt <- melt(
  combined_dt,
  id.vars = c("brain", "region", "file", "Object Id", "Analysis Region", "x_center", "x_rel"),
  measure.vars = c("ADCYAP1+", "RORB+", "PCP4+", "AQP4+"),
  variable.name = "celltype",
  value.name = "positive"
)

# Update colors to include Negative
full_cols <- c(
  "ADCYAP1+" = "#4cb048",
  "RORB+"    = "#994f9f",
  "PCP4+"    = "gold",
  "AQP4+"    = "pink",
  "Negative" = "grey90"
)


##### stacked bar plots for total cells by brain #######
## 1. Count unique total cells per sample (for our denominator)
#total_cells_dt <- plot_dt[, .(total_cells = uniqueN(`Object Id`)), by = .(brain, region)]
#
## 2. Count positive cells for each specific marker
#pos_counts <- plot_dt[positive == 1, .(N = .N), by = .(brain, region, celltype)]
#
## 3. Find cells that are negative for EVERYTHING
## We check each Object Id: if the sum of 'positive' is 0, it's a negative cell
#neg_counts <- plot_dt[, .(sum_pos = sum(positive)), by = .(brain, region, `Object Id`)]
#neg_counts <- neg_counts[sum_pos == 0, .(N = .N, celltype = "Negative"), by = .(brain, region)]
#
## 4. Combine them into one plotting table
#bar_plot_dt <- rbind(pos_counts, neg_counts, use.names = TRUE)
#
## 5. Calculate proportions
#bar_plot_dt <- merge(bar_plot_dt, total_cells_dt, by = c("brain", "region"))
#bar_plot_dt[, proportion := N / total_cells]
#
## Ensure 'Negative' is at the bottom of the stack for a cleaner look
#bar_plot_dt[, celltype := factor(celltype, levels = c("Negative", "AQP4+", "ADCYAP1+", "PCP4+", "RORB+"))]
#
#
#s_bar = ggplot(bar_plot_dt, aes(x = brain, y = proportion, fill = celltype)) +
#  # geom_col creates the bars; position = "stack" is the default
#  geom_col(color = "white", linewidth = 0.2) + 
#  facet_wrap(~region, ncol=1) +
#  scale_fill_manual(values = full_cols) +
#  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
#  labs(
#    x = "Sample (Donor)",
#    y = "Proportion of Total Cells",
#    fill = "Cell Category",
#    title = "Cell Type Composition by Region",
#    subtitle = "Comparing marker-positive vs. all-negative cells"
#  ) +
#  theme_classic(base_size = 14) +
#  theme(
#    strip.background = element_blank(),
#    strip.text = element_text(face = "bold", size = 15),
#    axis.text.x = element_text(angle = 45, hjust = 1),
#    legend.position = "right"
#  )
#  
#  #### stacked bar plots for total cells by xbins #######
#  # 1. Identify which 'Object Id' is negative for all markers
#  # We group by the cell's unique ID and check if it has zero positive markers
#  cell_status <- plot_dt[, .(is_any_pos = sum(positive) > 0), 
#                         by = .(region, x_bin, `Object Id`)]
#
#  # 2. Count "Negative" cells per bin/region
#  neg_bin_counts <- cell_status[is_any_pos == FALSE, .(N = .N), 
#                                 by = .(region, x_bin)]
#  neg_bin_counts[, celltype := "Negative"]
#
#  # 3. Count "Positive" cells per marker per bin/region
#  pos_bin_counts <- plot_dt[positive == 1, .(N = .N), 
#                             by = .(region, x_bin, celltype)]
#
#  # 4. Combine and calculate proportions per bin
#  bin_composition <- rbind(pos_bin_counts, neg_bin_counts, use.names = TRUE)
#
#  # Calculate the total cells in each bin to get the denominator
#  # This ensures each bar in the plot stacks exactly to 1.0 (100%)
#  bin_composition[, total_in_bin := sum(N), by = .(region, x_bin)]
#  bin_composition[, proportion := N / total_in_bin]
#
#  # 5. Set factor levels for a clean visual stack
#  bin_composition[, celltype := factor(celltype, levels = c("Negative", "AQP4+", "ADCYAP1+", "PCP4+", "RORB+"))]
#									   
#t_bar = ggplot(bin_composition, aes(x = x_bin, y = proportion, fill = celltype)) +
#    # width = 1 makes the bars touch, creating a continuous look
#    geom_col(width = 1) + 
#    facet_wrap(~region, ncol = 1) +
#    scale_fill_manual(values = full_cols) +
#    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
#    labs(
#      x = "Relative Depth (L1 -> WM)",
#      y = "Proportion of Cells",
#      fill = "Cell Type",
#      title = "Spatial Cell Type Composition across X-Bins",
#      subtitle = "Stacked proportions including Negative (all-marker 0) cells"
#    ) +
#    theme_classic(base_size = 14) +
#    theme(
#      axis.text.x = element_blank(), # Bins are 1-100, usually cleaner to hide labels
#      axis.ticks.x = element_blank(),
#      strip.background = element_blank(),
#      strip.text = element_text(face = "bold", size = 16),
#      panel.spacing = unit(1, "lines")
#    )
	
#### stacked bar plots for true total cells by xbins #######	
# 1. Determine the 'Status' for every single cell
#cell_status <- plot_dt[, .(
#  pos_count = sum(positive),
#  # If pos_count is 1, get the name of that marker; otherwise Label accordingly
#  category = ifelse(sum(positive) == 0, "Negative",
#             ifelse(sum(positive) > 1, "Multi-Positive", 
#                    as.character(celltype[positive == 1])))
#), by = .(region, x_bin, `Object Id`)]
#
# 2. Now count these mutually exclusive categories
#bin_composition <- cell_status[, .(N = .N), by = .(region, x_bin, category)]
#
## 3. Normalize by physical cell count
#bin_composition[, total_cells := sum(N), by = .(region, x_bin)]
#bin_composition[, proportion := N / total_cells]
#
#bin_composition[, category := factor(category, levels = c("Negative", "Multi-Positive", "AQP4+", "ADCYAP1+", "PCP4+", "RORB+"))]
## 4. Plot
#Tt_bar = ggplot(bin_composition, aes(x = x_bin, y = proportion, fill = category)) +
#  geom_col(width = 1) +
#  facet_wrap(~region, ncol = 1) +
#  scale_fill_manual(values = c(full_cols, "Multi-Positive" = "black")) + 
#  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
#  labs(
#    x = "Relative Depth (L1 -> WM)",
#    y = "Proportion of Cells",
#    fill = "Cell Type",
#    title = "True Spatial Cell Type Composition across X-Bins",
#    subtitle = "Stacked proportions including Negative (all-marker 0) cells"
#  ) +
#  theme_classic(base_size = 14) +
#  theme(  
#    axis.text.x = element_blank(), # Bins are 1-100ssource("daccVSdlpfc.R", chdir = TRUE)
#    axis.ticks.x = element_blank(),
#    strip.background = element_blank(),
#    strip.text = element_text(face = "bold", size = 16),
#    panel.spacing = unit(1, "lines")
#  )
  
##### density plots ########
# 1. Classify every unique cell based on its marker profile
cell_class_dt <- plot_dt[, .(
  pos_sum = sum(positive),
  # Find which marker is positive (only matters if pos_sum == 1)
  marker = as.character(celltype[positive == 1][1]) 
), by = .(brain, region, file, `Object Id`, x_rel)]

# 2. Assign the final mutually exclusive label
cell_class_dt[, category := "Negative"]
cell_class_dt[pos_sum == 1, category := marker]
cell_class_dt[pos_sum >= 2, category := "Multi-Positive"]

# 3. Get total object counts per region for weighting
region_totals <- cell_class_dt[, .(total_cells = .N), by = .(region)]
cell_class_dt <- merge(cell_class_dt, region_totals, by = "region")

# 4. Set factor levels for a logical plot order
cat_levels <- c("Negative", "Multi-Positive", "AQP4+", "ADCYAP1+", "PCP4+", "RORB+")
cell_class_dt[, category := factor(category, levels = cat_levels)]


#a_density = ggplot(cell_class_dt, aes(x = x_rel, fill = category, color = category, weight = 1/total_cells)) +
#  # Use a small alpha for Negative so it doesn't drown out the markers
#  geom_density(alpha = 0.4, linewidth = 0.6, position = "identity") +
#  facet_wrap(~region, ncol = 1) +
#  scale_fill_manual(values = c(full_cols, "Multi-Positive" = "black")) +
#  scale_color_manual(values = c(full_cols, "Multi-Positive" = "black")) +
#  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
#  labs(
#    x = "Relative Depth (0 = L1, 1 = WM)",
#    y = "Abundance Density",
#    title = "Spatial Density of Mutually Exclusive Cell Categories",
#    subtitle = "Area under curves sums to 1.0 (100% of cells)",
#    fill = "Classification",
#    color = "Classification"
#  ) +
#  theme_classic(base_size = 14) +
#  theme(strip.text = element_text(face = "bold"))
#  
#plot_dt_P <- plot_dt[positive == 1]
#
#p_density = ggplot(plot_dt_P, aes(x = x_rel, fill = celltype, color = celltype)) +
#  # alpha makes the overlapping areas visible
#  geom_density(alpha = 0.3, linewidth = 0.8) + 
#  scale_fill_manual(values = c(full_cols)) +
#  scale_color_manual(values = c(full_cols)) +
#  facet_wrap(~region, ncol = 1) +
#  labs(
#    x = "Relative Depth (0 = L1, 1 = WM)",
#    y = "Abundance Density",
#    title = "Spatial Density of Positive Cell types only",
#    subtitle = "Area under curves sums to 1.0 (100% of cells)",
#    fill = "Classification",
#    color = "Classification"
#  ) +
#  theme_classic(base_size = 14) +
#  theme(strip.text = element_text(face = "bold"))
#
#p_density

wewant = c("RORB+", "PCP4+", "ADCYAP1+")
cell_class_dt1 = cell_class_dt[cell_class_dt$category %in% wewant]

#aWe_density = ggplot(cell_class_dt1, aes(x = x_rel, fill = category, color = category, weight = 1/total_cells)) +
#  # Use a small alpha for Negative so it doesn't drown out the markers
#  geom_density(alpha = 0.4, linewidth = 0.6, position = "identity") +
#  facet_wrap(~region, ncol = 1) +
#  scale_fill_manual(values = c(full_cols)) +
#  scale_color_manual(values = c(full_cols)) +
#  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
#  labs(
#    x = "Relative Depth (0 = L1, 1 = WM)",
#    y = "Abundance Density",
#    title = "Spatial Density of Mutually Exclusive Cell Categories",
#    subtitle = "Area under curves sums to 1.0 (100% of cells)",
#    fill = "Classification",
#    color = "Classification"
#  ) +
#  theme_classic(base_size = 14) +
#  theme(strip.text = element_text(face = "bold"))

aWe_density_dACC = ggplot(cell_class_dt1[cell_class_dt1$region =="dACC", ], aes(y = x_rel, fill = category, color = category))+#, weight = 1/total_cells)) +
  # Use a small alpha for Negative so it doesn't drown out the markers
  geom_density(alpha = 0.4, linewidth = 0.6, position = "identity") +
  scale_fill_manual(values = c(full_cols)) +
  scale_color_manual(values = c(full_cols)) +
  #scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0), limits = c(1, 0)) +
  labs(
    y = "Relative Depth (0 = L1, 1 = WM)",
    x = "Abundance Density",
    #title = "Spatial Density of Mutually Exclusive Cell Categories",
    #subtitle = "Area under curves sums to 1.0 (100% of cells)",
    fill = "Classification",
    color = "Classification"
  ) +
  theme_classic(base_size = 14) +
  theme(strip.text = element_text(face = "bold"))

aWe_density_dlPFC = ggplot(cell_class_dt1[cell_class_dt1$region =="dlPFC", ], aes(y = x_rel, fill = category, color = category))+#, weight = 1/total_cells)) +
  # Use a small alpha for Negative so it doesn't drown out the markers
  geom_density(alpha = 0.4, linewidth = 0.6, position = "identity") +
  scale_fill_manual(values = c(full_cols)) +
  scale_color_manual(values = c(full_cols)) +
  #scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0), limits = c(1, 0)) +
  labs(
    y = "Relative Depth (0 = L1, 1 = WM)",
    x = "Abundance Density",
    #title = "Spatial Density of Mutually Exclusive Cell Categories",
    #subtitle = "Area under curves sums to 1.0 (100% of cells)",
    fill = "Classification",
    color = "Classification"
  ) +
  theme_classic(base_size = 14) +
  theme(strip.text = element_text(face = "bold"))


#plot_dt_P1<- plot_dt_P[plot_dt_P$celltype %in% wewant]
#pWe_density = ggplot(plot_dt_P1, aes(x = x_rel, fill = celltype, color = celltype)) +
#  # alpha makes the overlapping areas visible
#  geom_density(alpha = 0.3, linewidth = 0.8) + 
#  scale_fill_manual(values = c(full_cols)) +
#  scale_color_manual(values = c(full_cols)) +
#  facet_wrap(~region, ncol = 1) +
#  labs(
#    x = "Relative Depth (0 = L1, 1 = WM)",
#    y = "Abundance Density",
#    title = "Spatial Density of Positive Cell types only",
#    subtitle = "Area under curves sums to 1.0 (100% of cells)",
#    fill = "Classification",
#    color = "Classification"
#  ) +
#  theme_classic(base_size = 14) +
#  theme(strip.text = element_text(face = "bold"))
#
#pWe_density_dACC = ggplot(plot_dt_P1[plot_dt_P1$region == "dACC"], aes(y = x_rel, fill = celltype, color = celltype)) +
#  # alpha makes the overlapping areas visible
#  geom_density(alpha = 0.3, linewidth = 0.8) + 
#  scale_y_reverse(expand = c(0, 0), limits = c(1, 0)) +
#  scale_fill_manual(values = c(full_cols)) +
#  scale_color_manual(values = c(full_cols)) +
#  labs(
#    y = "Relative Depth (0 = L1, 1 = WM)",
#    x = "Abundance Density",
#    title = "Spatial Density of Positive Cell types only",
#    subtitle = "Area under curves sums to 1.0 (100% of cells)",
#    fill = "Classification",
#    color = "Classification"
#  ) +
#  theme_classic(base_size = 14) +
#  theme(strip.text = element_text(face = "bold"))
#
#pWe_density_dlPFC = ggplot(plot_dt_P1[plot_dt_P1$region == "dlPFC"], aes(y = x_rel, fill = celltype, color = celltype)) +
#  # alpha makes the overlapping areas visible
#  geom_density(alpha = 0.3, linewidth = 0.8) + 
#  scale_y_reverse(expand = c(0, 0), limits = c(1, 0)) +
#  scale_fill_manual(values = c(full_cols)) +
#  scale_color_manual(values = c(full_cols)) +
#  labs(
#    y = "Relative Depth (0 = L1, 1 = WM)",
#    x = "Abundance Density",
#    title = "Spatial Density of Positive Cell types only",
#    subtitle = "Area under curves sums to 1.0 (100% of cells)",
#    fill = "Classification",
#    color = "Classification"
#  ) +
#  theme_classic(base_size = 14) +
#  theme(strip.text = element_text(face = "bold"))
#
#  pWe_density_dACC <- pWe_density_dACC + theme(aspect.ratio = 2)
#  pWe_density_dlPFC <- pWe_density_dlPFC + theme(aspect.ratio = 2)
#
#  # 1. Define the output file path
#  pdf_out <- "/Users/madhavi.tippani/Desktop/HALO_dACC/HALO_Spatial_Analysis_Report.pdf"
#
#  # 2. Open the PDF device
#  # Using width=10 and height=8 is generally good for faceted horizontal plots
#  pdf(pdf_out, width = 10, height = 8)
#
#  # 3. Print the "Horizontal" plots (Stacked Bars and Standard Densities)
#  print(s_bar)
#  print(t_bar)
#  print(Tt_bar)
#  print(a_density)
#  print(p_density)
#  print(aWe_density)
#  print(pWe_density)
#  # 4. Print the "Vertical" plots
#  # Note: Since these are vertical, you might want a taller aspect ratio. 
#  # You can call pdf() again with new dimensions, or just print them here.
#  print(pWe_density_dACC)
#  print(pWe_density_dlPFC)
#
#  # 5. Close the device (CRITICAL: if you miss this, the PDF will be corrupted)
#  dev.off()
  
  
  aWe_density_dACC <- aWe_density_dACC + theme(aspect.ratio = 2)
  aWe_density_dlPFC <- aWe_density_dlPFC + theme(aspect.ratio = 2)

  # 1. Define the output file path
  pdf_out <- "/Users/madhavi.tippani/Desktop/HALO_dACC/HALO_Spatial_Analysis_density.pdf"

  # 2. Open the PDF device
  # Using width=10 and height=8 is generally good for faceted horizontal plots
  pdf(pdf_out, width = 10, height = 8)

  # 3. Print the "Horizontal" plots (Stacked Bars and Standard Densities)
  print(aWe_density_dACC)
  print(aWe_density_dlPFC)
 
  dev.off()
  
  # 1. Re-calculate centers for the whole dataset to include Y
  combined_dt[, y_center := (YMin + YMax) / 2]

  # 2. Subset for the specific brain
  brain_dt <- combined_dt[brain == "Br8667"]

  # 3. Get the "Positive" entries (Long format)
  # This naturally gives 2 rows for double-positive cells
  pos_spatial <- melt(
    brain_dt,
    id.vars = c("region", "x_center", "y_center", "Object Id"),
    measure.vars = c("ADCYAP1+", "RORB+", "PCP4+", "AQP4+"),
    variable.name = "celltype",
    value.name = "positive"
  )[positive == 1]

  # 4. Get the "Negative" entries (1 row per all-negative cell)
  # We find Object Ids that don't appear in the 'pos_spatial' set
  all_ids <- unique(brain_dt$`Object Id`)
  pos_ids <- unique(pos_spatial$`Object Id`)
  neg_ids <- setdiff(all_ids, pos_ids)

  neg_spatial <- brain_dt[`Object Id` %in% neg_ids, .(region, x_center, y_center, `Object Id`)]
  neg_spatial[, celltype := "Negative"]

  # 5. Combine for plotting
  spatial_plot_dt <- rbind(pos_spatial[, .(region, x_center, y_center, celltype)], 
                           neg_spatial[, .(region, x_center, y_center, celltype)])
						   
  # 1. Combine data: Negatives first (back layer), then Positives (front layer)
  spatial_plot_dt <- rbind(
    neg_spatial[, .(region, x_center, y_center, celltype)], 
    pos_spatial[, .(region, x_center, y_center, celltype)]
  )
             
  # 2. Set factor levels for the legend/colors
  spatial_plot_dt[, celltype := factor(celltype, levels = c("Negative", "AQP4+", "ADCYAP1+", "PCP4+", "RORB+"))]	
  spatial_plot_dt1 = spatial_plot_dt[spatial_plot_dt$celltype != "AQP4+", ]						   
  spatial_dACC <- ggplot(spatial_plot_dt1[spatial_plot_dt1$region == "dACC",], aes(y = x_center, x = y_center, color = celltype)) +
  # We use a very small alpha for the back layer (Negative) 
  # and a higher alpha for the front layer (Markers)
  geom_point(aes(alpha = celltype), size = 3, 
             position = position_jitter(width = 1.5, height = 1.5)) +
  scale_color_manual(values = full_cols) +
  #scale_y_reverse(expand = c(0, 0), limits = c(1, 0)) +
  # Set custom alpha: Negatives are faint (0.2), Markers are solid (0.9)
  scale_alpha_manual(values = c("Negative" = 0.7, "ADCYAP1+" = 0.9, "PCP4+" = 0.7, "RORB+" = 0.9)) +

  coord_fixed() + 
  theme_void() + 
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 16),
    panel.background = element_rect(fill = "black", color = NA)
  )

spatial_dACC

spatial_dlPFC <- ggplot(spatial_plot_dt1[spatial_plot_dt1$region == "dlPFC",], aes(y = x_center, x = y_center, color = celltype)) +
# We use a very small alpha for the back layer (Negative) 
# and a higher alpha for the front layer (Markers)
geom_point(aes(alpha = celltype), size = 3, 
           position = position_jitter(width = 1.5, height = 1.5)) +
scale_color_manual(values = full_cols) +
#scale_y_reverse(expand = c(0, 0), limits = c(1, 0)) +
# Set custom alpha: Negatives are faint (0.2), Markers are solid (0.9)
scale_alpha_manual(values = c("Negative" = 0.7, "ADCYAP1+" = 0.9, "PCP4+" = 0.5, "RORB+" = 0.9)) +

coord_fixed() + 
theme_void() + 
theme(
  legend.position = "none",
  strip.text = element_text(face = "bold", size = 16),
  panel.background = element_rect(fill = "black", color = NA)
)

spatial_dlPFC

spatial_dACC <- spatial_dACC + theme(aspect.ratio = 1.5)
spatial_dlPFC<- spatial_dlPFC + theme(aspect.ratio = 1.5)

# 1. Define the output file path
pdf_out <- "/Users/madhavi.tippani/Desktop/HALO_dACC/HALO_Spatial_Analysis_spatial.pdf"

# 2. Open the PDF device
# Using width=10 and height=8 is generally good for faceted horizontal plots
pdf(pdf_out, width = 12, height = 6)

# 3. Print the "Horizontal" plots (Stacked Bars and Standard Densities)
print(spatial_dACC)
print(spatial_dlPFC)

dev.off()



######## stats #########
# --- 4. T-STATS: PURE RORB+ vs PCP4+ (n=5 brains) ---
# Filter for only single-positive cells of interest

# 1. Isolate "Pure" Single-Positive Cells
#single_pos_dt <- plot_dt[, .(
#  pos_count = sum(positive),
#  phenotype = as.character(celltype[positive == 1][1])
#), by = .(brain, region, `Object Id`, x_rel)]
#
#pure_compare_dt <- single_pos_dt[pos_count == 1 & phenotype %in% c("RORB+", "PCP4+")]
#
## 2. Aggregate to Brain-Level Medians (n=5 per group)
#pure_brain_summary <- pure_compare_dt[, .(med_depth = median(x_rel)), 
#                                      by = .(brain, region, phenotype)]
#									  
pure_brain_summary <- cell_class_dt1[, .(med_depth = median(x_rel)), 
                                     by = .(brain, region, category)]
pure_brain_summary$phenotype = pure_brain_summary$category
pure_brain_summary[, phenotype := factor(phenotype, levels = c("RORB+", "PCP4+"))]

# 1. Run the Paired T-Test (Vector Method)
paired_stats <- pure_brain_summary[, {
  # Sort by brain to ensure the pairs match up correctly
  setorder(.SD, brain, phenotype)
  
  rorb_vals <- med_depth[phenotype == "RORB+"]
  pcp4_vals <- med_depth[phenotype == "PCP4+"]
  
  res <- t.test( pcp4_vals, rorb_vals, paired = TRUE)
  
  .(
    t_stat = res$statistic,
    p_value = res$p.value,
    mean_diff = res$estimate, # Mean of the differences
    conf_low = res$conf.int[1],
    conf_high = res$conf.int[2]
  )
}, by = region]

print("Paired T-test Results (n=5 pairs):")
print(paired_stats)

library(ggpubr)
tight_theme <- theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    # T, R, B, L margins - set to near zero or small values
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"),
    # Adjust aspect ratio if the boxes look too far apart
    aspect.ratio = 1.5 
  )
  
pure_paired_plot_dACC <- ggplot(pure_brain_summary[pure_brain_summary$region == "dACC"], aes(x = phenotype, y = med_depth, fill = phenotype)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA, width = 0.4) +
  
  # Connect the dots from the same brain
  geom_line(aes(group = brain), color = "grey70", alpha = 0.6) +
  geom_point(size = 3.5, shape = 21, color = "black") +
  
  #scale_y_reverse(limits = c(0.85, 0.15), expand = c(0, 0)) + 
 # scale_y_reverse(expand = expansion(mult = c(0.1, 0.2))) + 
  ylim(c(0.5, 0.9)) +
  scale_fill_manual(values = full_cols) +
  
  # Update to paired = TRUE
  stat_compare_means(method = "t.test", paired = TRUE, label = "p.format", label.y = 0.85, label.x =1.5, size = 5, color = "black") +
  
  labs(
    x = "Pure Phenotype",
    y = "Relative depth (L1 --> WM)",
    #title = "Laminar Separation: Paired Comparison",
    subtitle = "Lines connect RORB+ and PCP4+ medians from the same donor"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none", strip.text = element_text(face = "bold")) +tight_theme
  


pure_paired_plot_dlPFC <- ggplot(pure_brain_summary[pure_brain_summary$region == "dlPFC"], aes(x = phenotype, y = med_depth, fill = phenotype)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA, width = 0.4) +
  
  # Connect the dots from the same brain
  geom_line(aes(group = brain), color = "grey70", alpha = 0.6) +
  geom_point(size = 3.5, shape = 21, color = "black") +
  
  #scale_y_reverse(limits = c(0.85, 0.15), expand = c(0, 0)) + 
 # scale_y_reverse(expand = expansion(mult = c(0.1, 0.2))) + 
  ylim(c(0.5, 0.9)) +
  scale_fill_manual(values = full_cols) +
  
  # Update to paired = TRUE
  stat_compare_means(method = "t.test", paired = TRUE, label = "p.signif", label.y = 0.85, label.x =1.5, size = 10, color = "black") +
  
  labs(
    x = "Pure Phenotype",
    y = "Relative depth (L1 --> WM)",
    #title = "Laminar Separation: Paired Comparison",
    subtitle = "Lines connect RORB+ and PCP4+ medians from the same donor"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))+tight_theme
  
  
  pdf_out <- "/Users/madhavi.tippani/Desktop/HALO_dACC/HALO_Spatial_Analysis_stats.pdf"
  pdf(pdf_out, width = 8, height = 6)

  # 3. Print the "Horizontal" plots (Stacked Bars and Standard Densities)
  print(pure_paired_plot_dACC)
  print(pure_paired_plot_dlPFC)

  dev.off()


 