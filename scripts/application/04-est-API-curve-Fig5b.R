library(dplyr)
library(ggplot2)
# Load model outputs
out.summary <- readRDS("data/mdl-fit-outputs-summary.rds")

# Plot the estimated infection probability over time with credible intervals
est_s <- out.summary |>
  filter(grepl("^s_hat", variable)) |> 
  mutate(time = 1:100) |>
  mutate(year = 1917:2016) |> 
  rename(s_hat = mean) |> 
  ggplot() +
  # estimated curve
  geom_line(aes(x = year, y = s_hat, color = "Posterior mean"),
            linewidth = 1) +
  geom_ribbon(aes(ymin = q25, ymax = q75, x = year, alpha = "50%"), fill = "#2171B5") +
  geom_ribbon(aes(ymin = q10, ymax = q90, x = year, alpha = "80%"), fill = "#2171B5") +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5, x = year, alpha = "95%"), fill = "#2171B5") +
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
  ) + scale_x_continuous(
    breaks = c(1917, seq(1920, 2010, by = 10), 2016)
  ) +
  coord_cartesian(ylim = c(0, 0.12), clip = "off")  +
  labs(
    fill = "", color = "",
    x = "Time (year)",
    y = "Annual Probability of Infection",
    title = "" #"Estimation of infection probability over time"
  ) +
  theme_minimal() +                 
  theme(
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "darkgrey", linewidth = 0.5, linetype = "solid"),
    legend.position = c(0.9, 0.8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# ggsave("plots/Fig5b-est-api-curve-cis.pdf", plot = est_s, width = 10, height = 5)

