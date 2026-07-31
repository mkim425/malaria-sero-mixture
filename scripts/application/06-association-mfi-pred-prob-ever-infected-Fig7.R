library(dplyr)
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


p_mfi_prob <- data.with.prob |> 
  ggplot(aes(x = log(PvAMA1), y = prob_infection)) +
  geom_point(alpha = 0.4) +
  labs(
    title = "", #"Association between log(PvAMA1) and cumulative infection probability",
    x = "log(PvAMA1)",
    y = "Predicted infection probability"
  ) +
  theme_minimal() + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12)
  )

ggsave("plots/Fig7-mfi_prob.pdf",
       plot = p_mfi_prob, width = 10, height = 5)
