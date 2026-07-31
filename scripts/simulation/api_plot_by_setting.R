# function to plot s(t): annual probability of infection
api_plot_by_setting <- function(data1 = deterministic_df, 
                     data2 = spline_shat, 
                     data3 = coverage_df,
                     setting = "Setting 1"){
  ggplot() +
    # estimated curve
    geom_line(data = data2, aes(x = time, y = s_hat, group = sim_id),
              color = alpha("lightblue", 0.4), size = 0.7) +
    # deterministic curve
    geom_line(data = data1, aes(x = time, y = s_hat),
              color = "red", size = 0.7) +
    scale_x_continuous(
      breaks = seq(0, 100, by = 5)  # 5, 10, 15, 20
    ) +
    coord_cartesian(ylim = c(0, 0.105), clip = "off")  +
    labs(
      x = "",
      y = "",
      title = setting
    ) +
    theme_minimal() +                 
    theme(
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "darkgrey", linewidth = 0.5, linetype = "solid")
    ) +
    geom_text(data = data3,
              aes(x = time, y = 0.1, label = paste0(round(avg_coverage_95ci * 100,2), "%")),
              size = 2.75, angle=90,  color = "blue") + 
    annotate("text" , x = 110, y = 0.1, label = "95CI \n coverage", 
                 color = "blue", size = 2.75, fontface = "bold")
}
