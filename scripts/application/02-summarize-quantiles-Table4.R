library(dplyr)
library(sn)
library(posterior)

# Load model outputs
draws <- readRDS("data/mdl-fit-outputs-draws.rds")

# marginalized scale and shape parameters for the ever-infected group
draws$scale_mrg <- sqrt(draws[,"scale[2]"]^2 + draws[,"eta"]^2)
draws$shape_mrg <- (draws[,"scale[2]"] * draws[,"shape[2]"]) / sqrt(draws[,"scale[2]"]^2 + draws[,"eta"]^2 * (1 + draws[,"shape[2]"]^2))

# quantiles for the mixture component
q_values <- lapply(1:12000, function(i){
  q_comp1 <- sn::qsn(
    c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99),
    xi = draws[i,"location[1]"], 
    omega    = draws[i,"scale[1]"], 
    alpha    = draws[i,"shape[1]"]
  )
  
  q_comp2 <- sn::qsn(
    c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99),
    xi = draws[i,"location[2]"], 
    omega    = draws[i,"scale_mrg"], 
    alpha    = draws[i,"shape_mrg"]
  )
  
  data.frame(
    draw = i,
    q_comp1_1  = q_comp1[1],
    q_comp1_10 = q_comp1[2],
    q_comp1_25 = q_comp1[3],
    q_comp1_50 = q_comp1[4],
    q_comp1_75 = q_comp1[5],
    q_comp1_90 = q_comp1[6],
    q_comp1_99 = q_comp1[7],
    q_comp2_1  = q_comp2[1],
    q_comp2_10 = q_comp2[2],
    q_comp2_25 = q_comp2[3],
    q_comp2_50 = q_comp2[4],
    q_comp2_75 = q_comp2[5],
    q_comp2_90 = q_comp2[6],
    q_comp2_99 = q_comp2[7]
  )
}
) |> bind_rows()

# Summarize the quantiles
summary_statistics_2comps <- data.frame(
  variable = names(q_values)[-1],
  mean = sapply(q_values[-1], mean, na.rm = TRUE),
  l95  = sapply(q_values[-1], quantile, 0.025, na.rm = TRUE),
  u95  = sapply(q_values[-1], quantile, 0.975, na.rm = TRUE)
) |> 
  mutate(across(where(is.numeric), round, 2))
