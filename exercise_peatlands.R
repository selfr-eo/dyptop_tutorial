library(terra)
library(tidyterra)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(cowplot)
library(ggspatial)
library(dplyr)
library(tidyr)
library(purrr)
library(here)
library(khroma)
library(readr)

source(here::here("R/plot_discrete_cbar.R"))


rasta_p_over_pet <- rast(here("data/aridityindex_p_over_pet_zomeretal2022_v3_yr_1x1deg.nc"))

# resample to the same grid as files above
# Load the target raster (the reference grid)
target_raster <- rast(here("data/params_topmodel_M8_C12_filled_1x1deg.nc"))

# Regrid: Resample source raster to match the grid of the target raster
rasta_p_over_pet_regr <- resample(rasta_p_over_pet, target_raster, method = "bilinear")  # Options: "bilinear", "near", "cubic"

# create meaningful breaks, see https://geco-bern.github.io/les/ecohydrology.html#sec-budyko
breaks <- c(0, 0.03, 0.1, 0.2, 0.3, 0.4, 0.5, 0.65, 0.8, 1, 1.2, 1.6, 2, 2.5, 3, Inf)

# bin values to get a discrete color scale (personal preference)
rasta_p_over_pet_bin <- rasta_p_over_pet

values(rasta_p_over_pet_bin) <- cut(
  values(rasta_p_over_pet),
  breaks = breaks, 
  right = FALSE
)

# get coast outlines
layer_coast <- rnaturalearth::ne_coastline(
  scale = 110, 
  returnclass = "sf"
)

# get ocean layer
layer_ocean <- rnaturalearth::ne_download(
  scale = 110,
  type = "ocean",
  category = "physical",
  returnclass = "sf",
  destdir = here::here("data/")
)

# construct map
ggmap <- ggplot() +
  
  # aridity index raster
  tidyterra::geom_spatraster(
    data = rasta_p_over_pet_bin,
    show.legend = FALSE
  ) +
  
  # coastline
  geom_sf(
    data = layer_coast,
    colour = 'black',
    linewidth = 0.1
  ) +
  
  # ocean to mask zeros
  geom_sf(
    data = layer_ocean,
    color = NA,
    fill = "azure3"
  ) +
  
  # color palette from the khroma package
  scale_fill_roma(discrete = TRUE, name = "")  +
  coord_sf(
    ylim = c(-60, 85),
    expand = FALSE   # to draw map strictly bounded by the specified extent
  ) +
  xlab('') +
  ylab('') +
  theme_bw() +
  theme(axis.ticks.y.right = element_line(),
        axis.ticks.x.top = element_line(),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "white")
  )

gglegend <- plot_discrete_cbar(
  breaks = breaks,
  colors = c(khroma::color("roma")(length(breaks)-1)),
  legend_title = "",
  legend_direction = "vertical",
  width = 0.03,
  font_size = 3,
  expand_size_y = 0.5,
  spacing = "constant"
)

cowplot::plot_grid(ggmap, gglegend, ncol = 2, rel_widths = c(1, 0.10))


# soil carbon
# read water table depth outputs
# open raster
rasta_soilc <- rast(here("data/LPX-Bern_DYPTOP_vars_1990-2020_1x1deg_SOILC.nc"))

# extract data
df_soilc <- as.data.frame(rasta_soilc, xy = TRUE, na.rm = TRUE) |> 
  as_tibble()

# Data data is oddly organised in the NetCDF file. We therefore have to apply an
# additional step for re-organising it into a tidy format.
# Rename columns for clarity
colnames(df_soilc) <- c("lon", "lat", paste0("year_", 1990:2020))

# Convert from wide to long format (tidy)
df_soilc <- df_soilc |> 
  pivot_longer(cols = starts_with("year_"), names_to = "year", values_to = "soilc") |> 
  mutate(year = as.integer(gsub("year_", "", year)))  # Extract year

get_changerate <- function(df){
  linmod <- lm(soilc ~ year, data = df)  # Fit linear model: value ~ time
  coef(linmod)[2]                        # Extract slope (change rate)
}

get_mean <- function(df){
  mean(df$soilc)
}

df_soilc_nested <- df_soilc |> 
  group_by(lon, lat) |> 
  nest() |> 
  mutate(
    soilc_changerate = map_dbl(data, ~get_changerate(.)),
    soilc_mean = map_dbl(data, ~get_mean(.))
  ) |> 
  select(-data)   # to make object smaller again
# load global coastline data
world <- ne_coastline(scale = "small", returnclass = "sf")

# mean
soilc_mean <- ggplot() +
  
  # Add elevation layer
  geom_raster(
    data = df_soilc_nested, 
    aes(x = lon, y = lat, fill = soilc_mean*1e-3),
    show.legend = TRUE
  ) +
  scale_fill_batlowW(reverse = TRUE, name = expression(paste("Soil C \n(kgC m"^-2, ")"))) +  # Reverse the "lapaz" color scale
  theme_void() +
  theme(
    legend.position = "right", # Position the legend at the bottom of the plot
    legend.title = element_text(size = 10), # Adjust title font size
    legend.text = element_text(size = 8),    # Adjust legend text size
    panel.background = element_rect(fill = "grey70", color = NA)
  ) +
  coord_fixed() +
  geom_sf(data = world, fill = NA, color = "black", size = 0.05) +  # Continent outlines
  ylim(-55, 80) + 
  xlim(-180, 180)

# changerate
soilc_changerate <- ggplot() +
  
  # Add elevation layer
  geom_raster(
    data = df_soilc_nested, 
    aes(x = lon, y = lat, fill = soilc_changerate),
    show.legend = TRUE
  ) +
  scale_fill_vik(
    name = expression(paste("Soil C change \n(gC m"^-2, " yr"^-1, ")")),
    reverse = TRUE
  ) +
  theme_void() +
  theme(
    legend.position = "right", # Position the legend at the bottom of the plot
    legend.title = element_text(size = 10), # Adjust title font size
    legend.text = element_text(size = 8),    # Adjust legend text size
    panel.background = element_rect(fill = "grey70", color = NA)
  ) +
  coord_fixed() +
  geom_sf(data = world, fill = NA, color = "black", size = 0.05) +  # Continent outlines
  ylim(-55, 80) + 
  xlim(-180, 180)

plot_grid(soilc_mean, soilc_changerate, ncol = 1, labels = c("a", "b"))


# for creating labels
get_peat_crit_labels <- function(ai, soilc_changerate, soilc_mean){
  ifelse(
    ai > 1,
    ifelse(
      soilc_changerate > 10,
      "C accumulation ok",
      ifelse(
        soilc_mean > 50000,
        "C stock ok",
        "C stock too low")
    ),
    "too dry"
  )
}

# for creating booleans
get_peat_crit <- function(ai, soilc_changerate, soilc_mean){
  ifelse(
    ai > 1,
    ifelse(
      soilc_changerate > 10,
      TRUE,
      ifelse(
        soilc_mean > 50000,
        TRUE,
        FALSE)
    ),
    FALSE
  )
}


# merge data frame for potential peatland area fraction ...
df_combined <- df_soilc_nested |> 
  
  # ... and with P/PET
  left_join(
    as_tibble(rasta_p_over_pet_regr, xy = TRUE, na.rm = TRUE) |> 
      rename(lon = x, lat = y), 
    by = c("lon", "lat")
  )

df_combined <- df_combined |> 
  mutate(
    peat_crit_labels = get_peat_crit_labels(ai, soilc_changerate, soilc_mean),
    peat_crit = get_peat_crit(ai, soilc_changerate, soilc_mean)
  )



# plot
gg <- ggplot() +
  
  # Add elevation layer
  geom_raster(
    data = df_combined, 
    aes(x = lon, y = lat, fill = peat_crit_labels),
    show.legend = TRUE
  ) +
  theme_void() +
  theme(
    legend.position = "right", # Position the legend at the bottom of the plot
    legend.title = element_text(size = 10), # Adjust title font size
    legend.text = element_text(size = 8),    # Adjust legend text size
    panel.background = element_rect(fill = "grey70", color = NA)
  ) +
  coord_fixed() +
  geom_sf(data = world, fill = NA, color = "black", size = 0.05) +  # Continent outlines
  ylim(-55, 80) + 
  xlim(-180, 180) +
  scale_fill_manual(
    name = "",
    values = c(
      "too dry" = "wheat3",
      "C accumulation ok" = "dodgerblue4",
      "C stock ok" = "dodgerblue1",
      "C stock too low" = "wheat4"
    )
  )

gg

# potential peatland area
df_wtd_fflooded <- read_csv(here("data/df_wtd_fflooded.csv"))

take_third_highest <- function(vec){
  rev(sort(vec))[3]
}

df_peat <- df_wtd_fflooded |> 
  group_by(lon, lat) |> 
  summarise(fpeat_pot = take_third_highest(fflooded), .groups = "drop")

# potential peatland area fraction
gg <- ggplot() +
  
  # Add elevation layer
  geom_raster(
    data = df_peat, 
    aes(x = lon, y = lat, fill = fpeat_pot),
    show.legend = TRUE
  ) +
  scale_fill_batlowW(
    reverse = TRUE, 
    name = expression(paste(italic("f")[peat]^"pot"))
  ) +  # Reverse the "lapaz" color scale
  theme_void() +
  theme(
    legend.position = "right", # Position the legend at the bottom of the plot
    legend.title = element_text(size = 10), # Adjust title font size
    legend.text = element_text(size = 8),    # Adjust legend text size
    panel.background = element_rect(fill = "grey70", color = NA)
  ) +
  coord_fixed() +
  geom_sf(data = world, fill = NA, color = "black", size = 0.05) +  # Continent outlines
  ylim(-55, 80) + 
  xlim(-180, 180)

gg

# actual peatland

# merge data frame for potential peatland area fraction ...
df_combined2 <- df_combined |> 
  
  # ... with soil C outputs
  left_join(
    df_peat, 
    by = c("lon", "lat")
  )

get_fpeat_act <- function(fpeat_pot, peat_crit){
  ifelse(
    peat_crit,
    fpeat_pot,
    0
  )
}

df_combined2 <- df_combined2 |> 
  mutate(fpeat_act = get_fpeat_act(fpeat_pot, peat_crit))

# actual peatland area fraction
gg <- ggplot() +
  
  # Add elevation layer
  geom_raster(
    data = df_combined2, 
    aes(x = lon, y = lat, fill = fpeat_act),
    show.legend = TRUE
  ) +
  scale_fill_batlowW(
    reverse = TRUE, 
    name = expression(paste(italic("f")[peat]))
  ) +  # Reverse the "lapaz" color scale
  theme_void() +
  theme(
    legend.position = "right", # Position the legend at the bottom of the plot
    legend.title = element_text(size = 10), # Adjust title font size
    legend.text = element_text(size = 8),    # Adjust legend text size
    panel.background = element_rect(fill = "grey70", color = NA)
  ) +
  coord_fixed() +
  geom_sf(data = world, fill = NA, color = "black", size = 0.05) +  # Continent outlines
  ylim(-55, 80) + 
  xlim(-180, 180)

gg
