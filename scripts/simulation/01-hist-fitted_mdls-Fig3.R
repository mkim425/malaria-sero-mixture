# Plotting distributions by simulation setup

library(dplyr)
library(ggplot2)
library(bayesplot) # for ppc_dens_overlay
library(posterior)
library(sn) # for skew-normal distribution functions

# ----------------------------------------------------------------------------
# NOTE: 
# The mixture weight is individual-specific (depends on age via logtheta[n], n = 1,2,..., 1000), so
# for population-level plot, we use the average across individuals.
# For the ever-infected group, we use modified scale parameter and shape parameter 
# due to marginalization over the random effects in the model. 


# ----------------------------------------------------------------------------
## Setting 1: Two distributions are distinct 
# load the simulated data and the fitted model output (posterior draws of parameters) for setting 1
sim_df1 <- readRDS("data/simdata-setting1-sn-0001.rds")
draws_df1 <- readRDS("data/mdl-fit-simdata-setting1-0001-outputs-draws.rds")

set.seed(425)
# average theta across individuals for population-level mixture weight
logtheta_cols <- grep("logtheta", names(draws_df1), value = TRUE)
draws_df1$theta1 <- rowMeans(exp(draws_df1[, logtheta_cols]))  # avg P(never infected)
draws_df1$theta2 <- 1 - draws_df1$theta1  # avg P(ever infected)

# marginalized scale and shape parameters for the ever-infected group
draws_df1$scale_mrg <- sqrt(draws_df1[,"scale[2]"]^2 + draws_df1[,"eta"]^2)
draws_df1$shape_mrg <- (draws_df1[,"scale[2]"] * draws_df1[,"shape[2]"]) / sqrt(draws_df1[,"scale[2]"]^2 + draws_df1[,"eta"]^2 * (1 + draws_df1[,"shape[2]"]^2))

# evaluate densities on a grid for 500 draws
x_grid <- seq(min(log(sim_df1$PvAMA1)), max(log(sim_df1$PvAMA1)), length.out = 300)

rnd_idx <- sample(nrow(sim_df1), 500) 
curve_df1 <- lapply(rnd_idx, function(i) {
  data.frame(
    x      = x_grid,
    comp1  = dsn(x_grid, draws_df1[i,"location[1]"], draws_df1[i,"scale[1]"],  draws_df1[i,"shape[1]"])  * draws_df1[i,"theta1"],
    comp2  = dsn(x_grid, draws_df1[i,"location[2]"], draws_df1[i,"scale_mrg"], draws_df1[i,"shape_mrg"]) * draws_df1[i,"theta2"],
    draw   = i
  )
}) |> bind_rows()

p.setting1 <- ggplot() +
  geom_histogram(data = data.frame(y = log(sim_df1$PvAMA1)),
                 aes(x = y, y = after_stat(density)),
                 breaks = seq(min(log(sim_df1$PvAMA1)), max(log(sim_df1$PvAMA1)), length.out = 101),
                 fill = "grey80", color = "black") +
  geom_line(data = curve_df1, aes(x = x, y = comp1, group = draw, color = "Never infected"), alpha = 0.05) +
  geom_line(data = curve_df1, aes(x = x, y = comp2, group = draw, color = "Ever infected"), alpha = 0.05) +
  scale_color_manual(name = "", 
                     values = c("Never infected" = "blue", "Ever infected" = "red"),
                     breaks = c("Never infected", "Ever infected")) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = expression(log(y^{sim})), y = "Density", title = "(a) Setting 1") +
  theme_bw() + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold"),
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.text         = element_text(size = 10),
        legend.key.size     = unit(0.7, "cm")
        )


# ----------------------------------------------------------------------------
## Setting 2: Two distributions with some overlap
# load the simulated data and the fitted model output (posterior draws of parameters) for setting 2
sim_df2 <- readRDS("data/simdata-setting2-sn-0001.rds")
draws_df2 <- readRDS("data/mdl-fit-simdata-setting2-0001-outputs-draws.rds")

set.seed(425)

# average theta across individuals for population-level mixture weight
logtheta_cols <- grep("logtheta", names(draws_df2), value = TRUE)
draws_df2$theta1 <- rowMeans(exp(draws_df2[, logtheta_cols]))  # avg P(never infected)
draws_df2$theta2 <- 1 - draws_df2$theta1  # avg P(ever infected)

# marginalized scale and shape parameters for the ever-infected group
draws_df2$scale_mrg <- sqrt(draws_df2[,"scale[2]"]^2 + draws_df2[,"eta"]^2)
draws_df2$shape_mrg <- (draws_df2[,"scale[2]"] * draws_df2[,"shape[2]"]) / sqrt(draws_df2[,"scale[2]"]^2 + draws_df2[,"eta"]^2 * (1 + draws_df2[,"shape[2]"]^2))

# evaluate densities on a grid for 500 draws
x_grid <- seq(min(log(sim_df2$PvAMA1)), max(log(sim_df2$PvAMA1)), length.out = 300)

rnd_idx <- sample(nrow(sim_df2), 500) 
curve_df2 <- lapply(rnd_idx, function(i) {
  data.frame(
    x      = x_grid,
    comp1  = dsn(x_grid, draws_df2[i,"location[1]"], draws_df2[i,"scale[1]"],  draws_df2[i,"shape[1]"])  * draws_df2[i,"theta1"],
    comp2  = dsn(x_grid, draws_df2[i,"location[2]"], draws_df2[i,"scale_mrg"], draws_df2[i,"shape_mrg"]) * draws_df2[i,"theta2"],
    draw   = i
  )
}) |> bind_rows()

p.setting2 <- ggplot() +
  geom_histogram(data = data.frame(y = log(sim_df2$PvAMA1)),
                 aes(x = y, y = after_stat(density)),
                 breaks = seq(min(log(sim_df2$PvAMA1)), max(log(sim_df2$PvAMA1)), length.out = 101),
                 fill = "grey80", color = "black") +
  geom_line(data = curve_df2, aes(x = x, y = comp1, group = draw), color = "blue",  alpha = 0.05) +
  geom_line(data = curve_df2, aes(x = x, y = comp2, group = draw), color = "red",   alpha = 0.05) +
  labs(x = expression(log(y^{sim})), y = "Density", title = "(b) Setting 2") +
  theme_bw() + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold"))


# ----------------------------------------------------------------------------
## Setting 3: Two distributions with some overlap and individual-level random effects
# load the simulated data and the fitted model output (posterior draws of parameters) for setting 3
sim_df3 <- readRDS("data/simdata-setting3-sn-0001.rds")
draws_df3 <- readRDS("data/mdl-fit-simdata-setting3-0001-outputs-draws.rds")

set.seed(425)
# average theta across individuals for population-level mixture weight
logtheta_cols <- grep("logtheta", names(draws_df3), value = TRUE)
draws_df3$theta1 <- rowMeans(exp(draws_df3[, logtheta_cols]))  # avg P(never infected)
draws_df3$theta2 <- 1 - draws_df3$theta1  # avg P(ever infected)

# marginalized scale and shape parameters for the ever-infected group
draws_df3$scale_mrg <- sqrt(draws_df3[,"scale[2]"]^2 + draws_df3[,"eta"]^2)
draws_df3$shape_mrg <- (draws_df3[,"scale[2]"] * draws_df3[,"shape[2]"]) / sqrt(draws_df3[,"scale[2]"]^2 + draws_df3[,"eta"]^2 * (1 + draws_df3[,"shape[2]"]^2))

# evaluate densities on a grid for 500 draws
x_grid <- seq(min(log(sim_df3$PvAMA1)), max(log(sim_df3$PvAMA1)), length.out = 300)

rnd_idx <- sample(nrow(sim_df3), 500) 
curve_df3 <- lapply(rnd_idx, function(i) {
  data.frame(
    x      = x_grid,
    mix   = dsn(x_grid, draws_df3[i,"location[1]"], draws_df3[i,"scale[1]"],  draws_df3[i,"shape[1]"])  * draws_df3[i,"theta1"] +
      dsn(x_grid, draws_df3[i,"location[2]"], draws_df3[i,"scale_mrg"], draws_df3[i,"shape_mrg"]) * draws_df3[i,"theta2"],
    comp1  = dsn(x_grid, draws_df3[i,"location[1]"], draws_df3[i,"scale[1]"],  draws_df3[i,"shape[1]"])  * draws_df3[i,"theta1"],
    comp2  = dsn(x_grid, draws_df3[i,"location[2]"], draws_df3[i,"scale_mrg"], draws_df3[i,"shape_mrg"]) * draws_df3[i,"theta2"],
    draw   = i
  )
}) |> bind_rows()

p.setting3 <- ggplot() +
  geom_histogram(data = data.frame(y = log(sim_df3$PvAMA1)),
                 aes(x = y, y = after_stat(density)),
                 breaks = seq(min(log(sim_df3$PvAMA1)), max(log(sim_df3$PvAMA1)), length.out = 101),
                 fill = "grey80", color = "black") +
  geom_line(data = curve_df3, aes(x = x, y = comp1, group = draw), color = "blue",  alpha = 0.05) +
  geom_line(data = curve_df3, aes(x = x, y = comp2, group = draw), color = "red",   alpha = 0.05) +
  # geom_line(data = curve_df3, aes(x = x, y = mix,   group = draw), color = "black", alpha = 0.05) +
  labs(x = expression(log(y^{sim})), y = "Density", title = "(c) Setting 3") +
  theme_bw() + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold"))

# ----------------------------------------------------------------------------
library(ggpubr)
ggarrange(p.setting1, p.setting2, p.setting3,
          ncol = 3, nrow = 1)
# ggsave("plots/Fig3-simdata-hist-mdl-est-curves.pdf", width = 12, height = 5)









