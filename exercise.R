library(optimx)
library(terra)
library(tidyterra)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(cowplot)
library(ggspatial)
library(dplyr)
library(here)

source(here("R/plot_map.R"))



### 
# congo
file_path <- here("data/ga2_congo.nc")
rasta_cti <- rast(file_path)
bounding_box <- c(15,22,-3,3) # c(15,20,-2,3)  # max bb:  c(15,22,-3,3)
rasta_cti_cropped <- crop(rasta_cti, ext(bounding_box))

df_cti_congo <- as_tibble(rasta_cti_cropped, xy = TRUE, na.rm = TRUE) |>
  setNames(c("lon", "lat", "cti"))

gg1 <- df_cti_congo |>
  ggplot(aes(cti, ..density..)) +
  geom_density(color = "tomato") +
  labs(title = "Congo", x = "CTI", y = "density") +
  theme_classic() +
  xlim(0, NA)

gg2 <- df_cti_congo |>
  ggplot(aes(x = cti)) +
  stat_ecdf(geom = "step", color = "tomato") +  # ECDF plot
  labs(x = "CTI",
       y = "Cumulative density",
       title = " ") +
  theme_classic() +
  xlim(0, NA)

# switzerland
file_path <- here("data/ga2_switzerland.nc")
bounding_box <- c(xmin = 6, xmax = 10.5, ymin = 45.81799, ymax = 47.80838)
rasta_cti <- rast(file_path)
rasta_cti_cropped <- crop(rasta_cti, ext(bounding_box))

df_cti_switzerland <- as_tibble(rasta_cti_cropped, xy = TRUE, na.rm = TRUE) |>
  setNames(c("lon", "lat", "cti"))

gg3 <- df_cti_switzerland |>
  ggplot(aes(cti, ..density..)) +
  geom_density(color = "tomato") +
  labs(title = "Switzerland", x = "CTI", y = "density") +
  theme_classic() +
  xlim(0, NA)

gg4 <- df_cti_switzerland |>
  ggplot(aes(x = cti)) +
  stat_ecdf(geom = "step", color = "tomato") +  # ECDF plot
  labs(x = "CTI",
       y = "Cumulative density",
       title = " ") +
  theme_classic() +
  xlim(0, NA)

cowplot::plot_grid(gg1, gg2, gg3, gg4, ncol = 2)


#### p2 

m_par <- 8.0  # Stocker et al., 2014, Table 1
wtd <- -0.2   # m

calc_cti_crit <- function(vec_cti, wtd, m_par){
  mean(vec_cti) - m_par * wtd
}

# Congo
cti_crit_congo <- calc_cti_crit(df_cti_congo$cti, wtd, m_par)
cti_crit_congo
