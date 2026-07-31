library(dplyr)
library(ggplot2)
source("scripts/simulation/api_plot_by_setting.R")
# Load model fit results from setting 1
df_results_setting1 <- readRDS("data/mdl-fit-simdata-setting1-outputs-summary.rds")
df_results_setting2 <- readRDS("data/mdl-fit-simdata-setting2-outputs-summary.rds")
df_results_setting3 <- readRDS("data/mdl-fit-simdata-setting3-outputs-summary.rds")

# true s(t) curve
deterministic_df <- data.frame(time = 1:100) |>
  mutate(
    s_hat = 0.00055 + (0.055 - 0.00055) / (1 + 2 * exp(-0.1 * (70 - time)))
  )

true_s_curve <- deterministic_df |> 
  ggplot() +
  # deterministic curve
  geom_line(aes(x = time, y = s_hat),
            color = "red", linewidth = 1.5
  ) +
  scale_x_continuous(
    breaks = seq(0, 100, by = 5) # 5, 10, 15, 20
  ) +
  coord_cartesian(ylim = c(0, 0.105), clip = "off") +
  labs(
    x = "",
    y = "",
    title = "(a) True s(t) curve",
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "darkgrey", linewidth = 0.5, linetype = "solid")
  )

# Setting1-S hat (estimated annual probability of infection)
post_s_hat_summary_by_sim_setting1 <- df_results_setting1 |>
  filter(grepl("^s_hat", variable)) |>
  group_by(sim_id) |>
  mutate(time = 1:100) |>
  mutate(
    truth = 0.00055 + (0.055 - 0.00055) / (1 + 2 * exp(-0.1 * (70 - time))),
    coverage_90ci = ifelse(truth >= q5 & truth <= q95, 1, 0),
    coverage_95ci = ifelse(truth >= q2.5 & truth <= q97.5, 1, 0)
  ) |>
  ungroup()

coverage_s_hat_setting1 <- post_s_hat_summary_by_sim_setting1 |>
  group_by(time) |>
  summarise(
    avg_coverage_90ci = round(mean(coverage_90ci), 2),
    avg_coverage_95ci = round(mean(coverage_95ci), 2)
  )

s_hat_df_setting1 <- post_s_hat_summary_by_sim_setting1 |>
  select(sim_id, variable, mean, median, time)

coverage_df_sub_setting1 <- coverage_s_hat_setting1 |>
  filter(time %% 5 == 0)

p_setting1 <- api_plot_by_setting(deterministic_df,
                                  s_hat_df_setting1 |> rename(s_hat = mean),
                                  coverage_df_sub_setting1,
                                  setting = "(b) Setting 1"
)

# Setting2-S hat (estimated annual probability of infection)
post_s_hat_summary_by_sim_setting2 <- df_results_setting2 |>
  filter(grepl("^s_hat", variable)) |>
  group_by(sim_id) |>
  mutate(time = 1:100) |>
  mutate(
    truth = 0.00055 + (0.055 - 0.00055) / (1 + 2 * exp(-0.1 * (70 - time))),
    coverage_90ci = ifelse(truth >= q5 & truth <= q95, 1, 0),
    coverage_95ci = ifelse(truth >= q2.5 & truth <= q97.5, 1, 0)
  ) |>
  ungroup()

coverage_s_hat_setting2 <- post_s_hat_summary_by_sim_setting2 |>
  group_by(time) |>
  summarise(
    avg_coverage_90ci = round(mean(coverage_90ci), 2),
    avg_coverage_95ci = round(mean(coverage_95ci), 2)
  )

s_hat_df_setting2 <- post_s_hat_summary_by_sim_setting2 |>
  select(sim_id, variable, mean, median, time)

coverage_df_sub_setting2 <- coverage_s_hat_setting2 |>
  filter(time %% 5 == 0)

p_setting2 <- api_plot_by_setting(deterministic_df,
                                  s_hat_df_setting2 |> rename(s_hat = mean),
                                  coverage_df_sub_setting2,
                                  setting = "(c) Setting 2"
)


# Setting3-S hat (estimated annual probability of infection)
post_s_hat_summary_by_sim_setting3 <- df_results_setting3 |>
  filter(grepl("^s_hat", variable)) |>
  group_by(sim_id) |>
  mutate(time = 1:100) |>
  mutate(
    truth = 0.00055 + (0.055 - 0.00055) / (1 + 2 * exp(-0.1 * (70 - time))),
    coverage_90ci = ifelse(truth >= q5 & truth <= q95, 1, 0),
    coverage_95ci = ifelse(truth >= q2.5 & truth <= q97.5, 1, 0)
  ) |>
  ungroup()

coverage_s_hat_setting3 <- post_s_hat_summary_by_sim_setting3 |>
  group_by(time) |>
  summarise(
    avg_coverage_90ci = round(mean(coverage_90ci), 2),
    avg_coverage_95ci = round(mean(coverage_95ci), 2)
  )

s_hat_df_setting3 <- post_s_hat_summary_by_sim_setting3 |>
  select(sim_id, variable, mean, median, time)

coverage_df_sub_setting3 <- coverage_s_hat_setting3 |>
  filter(time %% 5 == 0)

p_setting3 <- api_plot_by_setting(deterministic_df,
                                  s_hat_df_setting3 |> rename(s_hat = mean),
                                  coverage_df_sub_setting3,
                                  setting = "(d) Setting 3"
)

# save
pdf("plots/Fig4-simdata-s_hat.pdf", width = 10, height = 8)
grid.arrange(
  arrangeGrob(true_s_curve,
              p_setting1,
              p_setting2,
              p_setting3,
              ncol = 2
  ),
  bottom = textGrob("Time (year)", gp = gpar(fontsize = 12)),
  left = textGrob("Annual Probability of Infection", rot = 90, gp = gpar(fontsize = 12))
)
dev.off()
