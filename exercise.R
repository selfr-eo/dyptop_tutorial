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

m_par <- 4.0  # Stocker et al., 2014, Table 1
wtd <- -0.2   # m

calc_cti_crit <- function(vec_cti, wtd, m_par){
  mean(vec_cti) - m_par * wtd
}

# Congo
cti_crit_congo <- calc_cti_crit(df_cti_congo$cti, wtd, m_par)
cti_crit_congo

# switzerland
cti_crit_switzerland <- calc_cti_crit(df_cti_switzerland$cti, wtd, m_par)
cti_crit_switzerland


calc_fflooded <- function(vec_cti, cti_crit){
  length(vec_cti[vec_cti > cti_crit])/length(vec_cti)
}

# Congo
fflooded_congo <- calc_fflooded(df_cti_congo$cti, cti_crit_congo)
fflooded_congo

# Switzerland
fflooded_switzerland <- calc_fflooded(df_cti_switzerland$cti, cti_crit_switzerland)
fflooded_switzerland

gg1 <- df_cti_congo |>
  ggplot(aes(x = cti)) +
  stat_ecdf(geom = "step", color = "tomato") +  # ECDF plot
  geom_vline(xintercept = cti_crit_congo, linetype = "dotted") +
  geom_hline(yintercept = 1 - fflooded_congo, linetype = "dotted") +
  labs(title = "Congo",
       x = "CTI",
       y = "Cumulative density") +
  theme_classic() +
  xlim(0, NA)

gg2 <- df_cti_switzerland |>
  ggplot(aes(x = cti)) +
  stat_ecdf(geom = "step", color = "tomato") +  # ECDF plot
  geom_vline(xintercept = cti_crit_switzerland, linetype = "dotted") +
  geom_hline(yintercept = 1 - fflooded_switzerland, linetype = "dotted") +
  labs(title = "Switzerland",
       x = "CTI",
       y = "Cumulative density") +
  theme_classic() +
  xlim(0, NA)

cowplot::plot_grid(gg1, gg2, ncol = 2)


# PLOT CTI MAPS

# congo
file_path <- here("data/ga2_congo.nc")
bounding_box <- c(15,22,-3,3)
map_congo <- plot_map(file_path, bounding_box, cti_crit = cti_crit_congo, show_legend = TRUE, show_inset = FALSE)
ggsave(here("book/images/map_congo_flooded.png"), plot = map_congo, width = 8, height = 8)

# switzerland
file_path <- here("data/ga2_switzerland.nc")
bounding_box <- c(xmin = 6, xmax = 10.5, ymin = 45.81799, ymax = 47.80838)
map_switzerland <- plot_map(file_path, bounding_box, cti_crit = cti_crit_switzerland, show_legend = TRUE, show_inset = FALSE)
ggsave(here("book/images/map_switzerland_flooded.png"), plot = map_switzerland, width = 8, height = 5)


# Flooding vs water table depth

calc_fflooded_wtd <- function(wtd, vec_cti, m_par){
  cti_crit <- calc_cti_crit(vec_cti, wtd, m_par)
  fflooded <- calc_fflooded(vec_cti, cti_crit)
  return(fflooded)
}

# Congo
df_wtd_congo <- tibble(wtd = seq(-1, 1, by = 0.01)) |>
  rowwise() |>
  mutate(fflooded = calc_fflooded_wtd(wtd, df_cti_congo$cti, m_par))

gg1 <- df_wtd_congo |>
  ggplot(aes(wtd, fflooded)) +
  geom_line() +
  theme_classic() +
  labs(title = "Congo", x = "Water table depth (m)", y = "Flooded fraction")

# Switzerland
df_wtd_switzerland <- tibble(wtd = seq(-1, 1, by = 0.01)) |>
  rowwise() |>
  mutate(fflooded = calc_fflooded_wtd(wtd, df_cti_switzerland$cti, m_par))

gg2 <- df_wtd_switzerland |>
  ggplot(aes(wtd, fflooded)) +
  geom_line() +
  theme_classic() +
  labs(title = "Switzerland", x = "Water table depth (m)", y = "Flooded fraction")

plot_grid(gg1, gg2, ncol = 2)

assymmetric_sigmoid <- function(x, par){
  ( 1 + par["v"] * exp(-par["k"] * (x - par["q"])))^(-1/par["v"])
}

# Generate example data
set.seed(1982)  # For reproducibility

# Define the loss function (SSE)
loss_function <- function(params, x, y) {
  par <- setNames(params, c("v", "k", "q"))  # Name parameters
  y_pred <- assymmetric_sigmoid(x, par)  # Compute predictions
  sum((y - y_pred)^2)  # Sum of squared errors
}

# Initial guesses for parameters
init_params <- c(v = 1, k = 0.5, q = 4)

# Congo
# Optimization using optimx
fit_congo <- optimx(
  par = init_params, 
  fn = loss_function, 
  x = df_wtd_congo$wtd, 
  y = df_wtd_congo$fflooded, 
  method = "BFGS"
)

# Extract best parameters
best_params <- coef(fit_congo)[1, ]
v_est_congo <- best_params["v"]
k_est_congo <- best_params["k"]
q_est_congo <- best_params["q"]

# Compute fitted values
df_congo <- data.frame(
  x = df_wtd_congo$wtd,
  y_obs = df_wtd_congo$fflooded,
  y_fit = assymmetric_sigmoid(df_wtd_congo$wtd, best_params)
)

# Plot results
gg1 <- ggplot(df_congo, aes(x = x)) +
  geom_point(aes(y = y_obs), color = "grey20", size = 2, alpha = 0.6) +  # Observations
  geom_line(aes(y = y_fit), color = "tomato", linewidth = 1) +  # Fitted curve
  labs(title = "Congo",
       subtitle = paste("Estimated v =", round(v_est_congo, 2), 
                        "k =", round(k_est_congo, 2),
                        "q =", round(q_est_congo, 2)),
       x = "Water table depth (m)", y = "Flooded fraction") +
  theme_classic()

# Switzerland
# Optimization using optimx
fit_switzerland <- optimx(
  par = init_params, 
  fn = loss_function, 
  x = df_wtd_switzerland$wtd, 
  y = df_wtd_switzerland$fflooded, 
  method = "BFGS"
)

# Extract best parameters
best_params <- coef(fit_switzerland)[1, ]
v_est_switzerland <- best_params["v"]
k_est_switzerland <- best_params["k"]
q_est_switzerland <- best_params["q"]

# Compute fitted values
df_switzerland <- data.frame(
  x = df_wtd_switzerland$wtd,
  y_obs = df_wtd_switzerland$fflooded,
  y_fit = assymmetric_sigmoid(df_wtd_switzerland$wtd, best_params)
)

# Plot results
gg2 <- ggplot(df_switzerland, aes(x = x)) +
  geom_point(aes(y = y_obs), color = "grey20", size = 2, alpha = 0.6) +  # Observations
  geom_line(aes(y = y_fit), color = "tomato", linewidth = 1) +  # Fitted curve
  labs(title = "Switzerland",
       subtitle = paste("Estimated v =", round(v_est_switzerland, 2), 
                        "k =", round(k_est_switzerland, 2),
                        "q =", round(q_est_switzerland, 2)),
       x = "Water table depth (m)", y = "Flooded fraction") +
  theme_classic()

plot_grid(gg1, gg2, ncol = 2)


