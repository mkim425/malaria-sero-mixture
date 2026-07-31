library(dplyr)
library(sn)
library(posterior)
library(ggplot2)
library(ggpubr)

# Load the Lao held-out inference data 
data <- readRDS("data/test_serology_sample.rds")
data$AGE_GP <- cut(
  data$AGE,
  breaks = c(seq(0, 70, by = 10), Inf),
  right = FALSE,
  labels = c(
    "[0,10)", "[10,20)", "[20,30)", "[30,40)",
    "[40,50)", "[50,60)", "[60,70)", "70+"
  )
)
# Load model outputs
out.draws <- readRDS("data/mdl-fit-outputs-draws.rds")

draws_y_tilde <- out.draws |>
  subset_draws(variable = c("y_tilde")) |> 
  as_draws_matrix() |> 
  as.data.frame()

# density overlay on the log scale using posterior predictive draws
set.seed(425)
draws_y_tilde_mat <- as.matrix(draws_y_tilde)
ppc_y_tilde_overlay <- ppc_dens_overlay(
  y = log(data$PvAMA1),
  yrep = draws_y_tilde_mat[sample(
    seq_len(nrow(draws_y_tilde_mat)),
    size = min(200, nrow(draws_y_tilde_mat))
  ), ]
) +
  coord_cartesian(xlim = c(-0.2, 12)) +
  labs(
    x = "log(PvAMA1)",
    y = "Density",
    title = ""
  ) +
  theme_bw() + 
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 22)
  )


ggsave("plots/Fig8a-ppc_y_tilde_overlay.pdf",
       plot =ppc_y_tilde_overlay, width = 10, height = 5)

# grouped PPC density overlays by age group
ppc_y_tilde_overlay_by_agegp <- lapply(levels(data$AGE_GP), function(age_gp) {
  obs_idx <- data$AGE_GP == age_gp
  y_obs_gp <- log(data$PvAMA1[obs_idx])
  yrep_gp <- draws_y_tilde_mat[, obs_idx, drop = FALSE]
  yrep_gp <- yrep_gp[sample(
    seq_len(nrow(yrep_gp)),
    size = min(200, nrow(yrep_gp))
  ), , drop = FALSE]
  
  ppc_dens_overlay(y = y_obs_gp, yrep = yrep_gp) +
    coord_cartesian(xlim = c(-0.2, 12)) +
    labs(
      x = "",
      y = "",
      title = paste("Age group", age_gp)
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(
        size = 10,      
        hjust = 0.5     
      ),
      legend.text = element_text(size = 12)
    )
})

ppc_y_tilde_overlay_by_agegp_combined <- ggarrange(
  plotlist = ppc_y_tilde_overlay_by_agegp,
  ncol = 2, nrow = 4,
  common.legend = TRUE,
  legend = "bottom"
)

ppc_y_tilde_overlay_by_agegp_final  <- annotate_figure(
  ppc_y_tilde_overlay_by_agegp_combined,
  left = text_grob("Density", rot = 90),
  bottom = text_grob("log(PvAMA1)")
)

ggsave("plots/Fig8b-ppc_y_tilde_overlay_by_agegp.pdf",
       plot =ppc_y_tilde_overlay_by_agegp_final, width = 8, height = 8)
