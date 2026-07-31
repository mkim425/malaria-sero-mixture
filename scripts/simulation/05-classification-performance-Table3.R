library(dplyr)

classification_metrics <- readRDS("data/mdl-fit-simdata-outputs-classification.rds")

res <- classification_metrics |>
  group_by(Setting) |> 
  summarise(
    mean_AUC = mean(AUC),
    sd_AUC = sd(AUC),
    mean_Sensitivity = mean(Sensitivity),
    sd_Sensitivity = sd(Sensitivity),
    mean_Specificity = mean(Specificity),
    sd_Specificity = sd(Specificity),
    mean_Accuracy = mean(Accuracy),
    sd_Accuracy = sd(Accuracy),
    mean_Brier = mean(Brier),
    sd_Brier = sd(Brier)
  ) |> 
  mutate(
    across(where(is.numeric), \(x) round(x, 3))
  )

