library(dplyr)
library(posterior)
library(ggplot2)


# Load the Lao held-out inference data 
data <- readRDS("data/test_serology_sample.rds")
# add an age group variable to the data
data$AGE_GP <- cut(
  data$AGE,
  breaks = seq(0, 100, by = 10),
  right = FALSE
)

# Load model outputs
out.draws <- readRDS("data/mdl-fit-outputs-draws.rds")

# Extract the draws for the probability of no infection
draws <- out.draws |>
  subset_draws(variable = c("prob_noinfection")) |> 
  as_draws_matrix() |> 
  as.data.frame()

prob_noinfection_draws <- draws[, grep("prob_noinfection", colnames(draws))]

# predicted probability of ever-infected = 1- prob_noinfection
prob_inf_mean  <- 1 - colMeans(prob_noinfection_draws) 

data.with.prob <- data.frame(
  prob_infection = prob_inf_mean
) |> 
  bind_cols(data) |> 
  as_tibble() |> 
  select(AGE, AGE_GP, PvAMA1, prob_infection)

# Figure 6(a)(b)----------------------------------------------------------------
p1 <- data.with.prob |> 
  ggplot() +
  geom_histogram(aes(x = prob_infection, y = after_stat(count / sum(count))),
                 breaks = seq(min(data.with.prob$prob_infection), 
                              max(data.with.prob$prob_infection), 
                              length.out = 50),
                 fill = "grey80", color = "black") +
  labs(
    x = "", # "Predicted probability of ever-infected",
    y = "" ,#"Density",
    title = "(a) Full range" #Distribution of predicted infection probability"
  ) +
  theme_minimal() +                 
  theme(
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "darkgrey", linewidth = 0.5, linetype = "solid"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 15)
  )

p2 <- data.with.prob |> 
  filter(prob_infection > 0.1, prob_infection < 0.9)  |> 
  ggplot() +
  geom_histogram(aes(x = prob_infection, y = after_stat(count / sum(count))), 
                 bins = 40,
                 fill = "grey80", color = "black") +
  labs(
    x = "", # "Predicted probability of ever-infected",
    y = "" ,#"Density",
    title = "(b) Zoomed region between 0.1 and 0.9" #"Distribution of predicted infection probability between 0.1 and 0.9 "
  ) +
  theme_minimal() +                 
  theme(
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "darkgrey", linewidth = 0.5, linetype = "solid"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 15)
  )

hist_prob_infection_prop <- grid.arrange(
  arrangeGrob(p1, p2, ncol = 2),
  top = textGrob("",
                 gp = gpar(fontsize = 14, fontface = "bold")
  ),
  bottom = textGrob("Predicted probability of ever-infected", gp = gpar(fontsize = 14)),
  left = textGrob("Proportion", rot = 90, gp = gpar(fontsize = 14))
)


# ggsave("plots/Fig6ab-hist_prob_infection_prop.pdf",
#        plot = hist_prob_infection_prop, width = 10, height = 5)


# Figure 6(c) ------------------------------------------------------------------
probinfected_agegp_prop <- data.with.prob |> 
  mutate(AGE_GP = ifelse(AGE_GP == "[90,100)", "70+", as.character(AGE_GP)),
         AGE_GP = ifelse(AGE_GP == "[80,90)", "70+", as.character(AGE_GP)),
         AGE_GP = ifelse(AGE_GP == "[70,80)", "70+", as.character(AGE_GP))) |>
  ggplot(aes(x = prob_infection)) +
  geom_histogram(aes(y = after_stat(count / sum(count))), 
                 binwidth = 0.05, fill = "grey80", color = "black") +
  facet_wrap(~AGE_GP, scales = "free_y") +
  labs(x = "Predicted probability of ever-infected", y = "Proportion",
       title = "(c) Distribution by age group" 
  ) +
  theme_bw() 

# ggsave("plots/Fig6c-probinfected_agegp.pdf", 
#        plot =probinfected_agegp_prop, width = 8, height = 6)



