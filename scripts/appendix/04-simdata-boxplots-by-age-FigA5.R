# Plotting box plots by simulation setup

library(dplyr)
library(ggplot2)

## Setting 1: Two distributions are distinct 
# load the simulated data for setting 1
sim_df1 <- readRDS("data/simdata-setting1-sn-0001.rds")
sim_df1$AGE_GP <- cut(
  sim_df1$age,
  breaks = seq(0, 100, by = 10),
  right = FALSE
)

bp1 <- sim_df1 |> 
  ggplot() +
  geom_boxplot(aes(x = AGE_GP, y = log(PvAMA1))) +
  labs(x = "Age group", y = expression(log(y^{sim})),
       title = "(a) Setting 1" 
  ) +
  theme_bw() + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# ----------------------------------------------------------------------------
## Setting 2: Two distributions with some overlap
# load the simulated data for setting 2
sim_df2 <- readRDS("data/simdata-setting2-sn-0001.rds")
sim_df2$AGE_GP <- cut(
  sim_df2$age,
  breaks = seq(0, 100, by = 10),
  right = FALSE
)

bp2 <- sim_df2 |> 
  ggplot() +
  geom_boxplot(aes(x = AGE_GP, y = log(PvAMA1))) +
  labs(x = "Age group", 
       y = expression(log(y^{sim})),
       title = "(b) Setting 2" 
  ) +
  theme_bw() + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# ----------------------------------------------------------------------------
## Setting 3: Two distributions with some overlap and individual-level random effects
# load the simulated data for setting 3
sim_df3 <- readRDS("data/simdata-setting3-sn-0001.rds")
sim_df3$AGE_GP <- cut(
  sim_df3$age,
  breaks = seq(0, 100, by = 10),
  right = FALSE
)

bp3 <- sim_df3 |> 
  ggplot() +
  geom_boxplot(aes(x = AGE_GP, y = log(PvAMA1))) +
  labs(x = "Age group", 
       y = expression(log(y^{sim})),
       title = "(c) Setting 3" 
  ) +
  theme_bw() + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())

library(ggpubr)
bp_all <- ggarrange(bp1, bp2, bp3,
          ncol = 1, nrow = 3, align = "v")
annotate_figure(bp_all, 
                bottom = textGrob("Age group", gp = gpar(cex = 1.2)))
# ggsave("plots/FigA5-simdata-boxplots-by-age.pdf", width = 6, height = 9)