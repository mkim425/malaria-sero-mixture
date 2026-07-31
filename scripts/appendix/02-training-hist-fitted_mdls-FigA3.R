library(dplyr)
library(tidyr)
library(ggplot2)
library(sn)
library(gridExtra)
library(grid)
# --- Lao model development data
data <- readRDS("data/training_serology_sample_clean.rds")

# --- model fit outputs
out.summary <- readRDS("data/mdl-fit-training-outputs-summary.rds")
out.draws <- readRDS("data/mdl-fit-training-outputs-draws-reduced.rds")
draws <- out.draws |> as.data.frame()

# set seed for reproducibility
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
                 breaks = seq(min(log(data$PvAMA1)), max(log(data$PvAMA1)), 
                              length.out = 101),
                 fill = "grey80", color = "black") +
  geom_line(data = curve, aes(x = x, y = comp1, group = draw, color = "Never infected"), alpha = 0.05) +
  geom_line(data = curve, aes(x = x, y = comp2, group = draw, color = "Ever infected"), alpha = 0.05) +
  scale_color_manual(name = "", 
                     values = c("Never infected" = "blue", "Ever infected" = "red"),
                     breaks = c("Never infected", "Ever infected")) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = expression(log(PvAMA1)), y = "Density", title = "(a)") +
  theme_bw() + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold"),
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.text         = element_text(size = 10),
        legend.key.size     = unit(0.7, "cm")
  )

est_s <- out.summary |>
  filter(grepl("^s_hat", variable)) |> 
  mutate(time = 1:100) |>
  # select(variable, mean, median, time) |> 
  rename(s_hat = mean) |> 
  ggplot() +
  # estimated curve
  geom_line(aes(x = time, y = s_hat, color = "Posterior mean"),
            linewidth = 1) +
  geom_ribbon(aes(ymin = q25, ymax = q75, x = time, alpha = "50%"), fill = "#2171B5") +
  geom_ribbon(aes(ymin = q10, ymax = q90, x = time, alpha = "80%"), fill = "#2171B5") +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5, x = time, alpha = "95%"), fill = "#2171B5") +
  scale_color_manual(
    values = c(
      "Posterior mean" = "blue"
    )
  ) + 
  scale_alpha_manual(
    name = "Credible interval",
    values = c(
      "95%" = 0.2,
      "80%" = 0.35,
      "50%" = 0.5
    )
  ) +
  scale_x_continuous(
    breaks = seq(0, 100, by = 5)  # 5, 10, 15, 20
  ) +
  coord_cartesian(ylim = c(0, 0.12), clip = "off")  +
  labs(
    fill = "", color = "",
    x = "Time (year)",
    y = "Annual Probability of Infection",
    title = "(b)"
  ) +
  theme_minimal() +                 
  theme(
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "darkgrey", linewidth = 0.5, linetype = "solid"),
    legend.position = c(0.8, 0.8)
  )

# combine plots
training_mdl_fit_result <- grid.arrange(
  arrangeGrob(hist_est_curve, est_s, ncol = 2, widths = c(0.4, 0.6)),
  top = textGrob("", #"Using model development data",
                 gp = gpar(fontsize = 14, fontface = "bold")
  )
)

# save
ggsave("plots/FigA3-training-hist-mdl-est-curve.pdf",
       plot = training_mdl_fit_result, width = 10, height = 5)
