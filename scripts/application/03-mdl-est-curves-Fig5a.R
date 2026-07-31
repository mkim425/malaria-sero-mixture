library(dplyr)
library(sn)
library(posterior)
library(ggplot2)

# Load the Lao held-out inference data 
data <- readRDS("data/test_serology_sample.rds")

# Load model outputs
out.draws <- readRDS("data/mdl-fit-outputs-draws.rds")

draws <- out.draws |>
  subset_draws(variable = c("location", "scale", "shape", "eta", "logtheta")) |> 
  as_draws_matrix() |> 
  as.data.frame()

set.seed(425)
# average theta across individuals for population-level mixture weight
logtheta_cols <- grep("logtheta", names(draws), value = TRUE)
draws$theta1 <- rowMeans(exp(draws[, logtheta_cols]))  # avg P(never infected)
draws$theta2 <- 1 - draws$theta1  # avg P(ever infected)

# marginalized scale and shape parameters for the ever-infected group
draws$scale_mrg <- sqrt(draws[,"scale[2]"]^2 + draws[,"eta"]^2)
draws$shape_mrg <- (draws[,"scale[2]"] * draws[,"shape[2]"]) / sqrt(draws[,"scale[2]"]^2 + draws[,"eta"]^2 * (1 + draws[,"shape[2]"]^2))

# evaluate densities on a grid for 500 draws
x_grid <- seq(min(log(data$PvAMA1)), max(log(data$PvAMA1)), length.out = 300)

rnd_idx <- sample(nrow(data), 500) 
curve <- lapply(rnd_idx, function(i) {
  data.frame(
    x      = x_grid,
    comp1  = dsn(x_grid, draws[i,"location[1]"], draws[i,"scale[1]"],  draws[i,"shape[1]"])  * draws[i,"theta1"],
    comp2  = dsn(x_grid, draws[i,"location[2]"], draws[i,"scale_mrg"], draws[i,"shape_mrg"]) * draws[i,"theta2"],
    draw   = i
  )
}) |> bind_rows()

hist_est_curve <- ggplot() +
  geom_histogram(data = data.frame(y = log(data$PvAMA1)),
                 aes(x = y, y = after_stat(density)),
                 breaks = seq(min(log(data$PvAMA1)), max(log(data$PvAMA1)), length.out = 101),
                 fill = "grey80", color = "black") +
  geom_line(data = curve, aes(x = x, y = comp1, group = draw, color = "Never infected"), alpha = 0.05) +
  geom_line(data = curve, aes(x = x, y = comp2, group = draw, color = "Ever infected"), alpha = 0.05) +
  scale_color_manual(name = "", 
                     values = c("Never infected" = "blue", "Ever infected" = "red"),
                     breaks = c("Never infected", "Ever infected")) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = expression(log(PvAMA1)), y = "Density", title = "") +
  theme_bw() + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold"),
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.text         = element_text(size = 10),
        legend.key.size     = unit(0.7, "cm")
  )

# ggsave("plots/Fig5a-hist-mdl-est-curve.pdf",
#         plot = hist_est_curve, width = 10, height = 5)
