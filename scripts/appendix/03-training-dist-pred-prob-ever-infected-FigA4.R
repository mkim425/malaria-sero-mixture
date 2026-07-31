library(dplyr)
library(posterior)
library(ggplot2)
# --- Lao model development data
data <- readRDS("data/training_serology_sample_clean.rds")
# add an age group variable to the data
data$AGE_GP <- cut(
  data$AGE,
  breaks = seq(0, 100, by = 10),
  right = FALSE
)
# --- model fit outputs
out.draws <- readRDS("data/mdl-fit-training-outputs-draws-prob_noinfection.rds") 

# predicted probability of ever-infected = 1- prob_noinfection
prob_inf_mean  <- 1 - colMeans(out.draws) 

data.with.prob <- data.frame(
  prob_infection = prob_inf_mean
) |> 
  bind_cols(data) |> 
  as_tibble() |> 
  select(AGE, AGE_GP, PvAMA1, prob_infection)

training_probinfected_agegp_prop <- data.with.prob |> 
  mutate(AGE_GP = ifelse(AGE_GP == "[90,100)", "70+", as.character(AGE_GP)),
         AGE_GP = ifelse(AGE_GP == "[80,90)", "70+", as.character(AGE_GP)),
         AGE_GP = ifelse(AGE_GP == "[70,80)", "70+", as.character(AGE_GP))) |>
  ggplot(aes(x = prob_infection)) +
  geom_histogram(aes(y = after_stat(count / sum(count))), 
                 binwidth = 0.05, fill = "grey80", color = "black") +
  facet_wrap(~AGE_GP, scales = "free_y") +
  labs(x = "Predicted probability of ever-infected", y = "Proportion",
       title = "" 
  ) +
  theme_bw() 

ggsave("plots/FigA4-training_probinfected_agegp.pdf",
       plot = training_probinfected_agegp_prop, width = 8, height = 6)
