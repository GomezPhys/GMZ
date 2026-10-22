# 1. Load packages
library(readxl)
library(dplyr)
library(tidyr)
library(rstatix)
library(ggplot2)


Peak <- read_excel("~/GXT_TTF.xlsx", sheet = "GXT_Peak")
View(Peak)

## Order conditions
Peak$Condition <- ordered(Peak$Condition,
                        levels = c("ECC", "CON"))

# 4. Define physiological variables
variables <- c(
  "VO2",
  "HR",
  "RPE",
  "Lactate",
  "PO",
  "S_BP",
  "D_BP"
)

# 5. Reshape data temporarily to wide format
# Each participant must have one observation per condition

Peak_wide <- Peak %>%
  select(ID, Condition, all_of(variables)) %>%
  pivot_wider(
    id_cols = ID,
    names_from = Condition,
    values_from = all_of(variables),
    names_glue = "{.value}_{Condition}"
  )

# 6. Calculate paired differences (CON - ECC)

for (var in variables) {
  Peak_wide[[paste0(var, "_diff")]] <-
    Peak_wide[[paste0(var, "_CON")]] -
    Peak_wide[[paste0(var, "_ECC")]]
}

# 7. Run Shapiro-Wilk tests for all variables

normality_results <- lapply(variables, function(var) {
  
  differences <- Peak_wide[[paste0(var, "_diff")]]
  differences <- differences[is.finite(differences)]
  
  if (length(differences) < 3 ||
      length(differences) > 5000 ||
      length(unique(differences)) < 2) {
    return(data.frame(
      Variable = var,
      N = length(differences),
      W = NA_real_,
      P_value = NA_real_,
      Interpretation = "Test not applicable"
    ))
  }
  
  test <- shapiro.test(differences)
  
  data.frame(
    Variable = var,
    N = length(differences),
    W = round(unname(test$statistic), 4),
    P_value = round(test$p.value, 4),
    Interpretation = ifelse(
      test$p.value > 0.05,
      "No significant deviation",
      "Significant deviation"
    )
  )
})

# 8. Combine results into one table

normality_results <- bind_rows(normality_results)

# Display results
print(normality_results)
View(normality_results)


# ============================================================
# STEP 2: Q-Q PLOTS
# ============================================================

# Convert paired differences to long format
Peak_differences <- Peak_wide %>%
  select(ID, ends_with("_diff")) %>%
  pivot_longer(
    cols = ends_with("_diff"),
    names_to = "Variable",
    values_to = "Difference"
  ) %>%
  mutate(
    Variable = sub("_diff$", "", Variable)
  ) %>%
  filter(is.finite(Difference))

# Generate Q-Q plots for all variables
ggplot(Peak_differences, aes(sample = Difference)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  facet_wrap(~ Variable, scales = "free") +
  theme_classic() +
  labs(
    title = "GXT Peak: Normality of Paired Differences",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  )

# Create cleaned VO2 differences
VO2_clean <- Peak_wide %>%
  filter(ID != "GMZ_11") %>%
  pull(VO2_diff)

# Repeat Shapiro-Wilk
shapiro.test(VO2_clean)

# Repeat Q-Q plot
qqnorm(VO2_clean)
qqline(VO2_clean, col = "red", lwd = 2)


####LMM#####
library(dplyr)
library(lmerTest)
library(emmeans)
library(ggplot2)

VO2_data <- Peak %>%
  filter(Condition %in% c("ECC", "CON")) %>%
  filter(!is.na(VO2)) %>%
  mutate(
    ID = factor(ID),
    Condition = factor(Condition, levels = c("ECC", "CON"))
  )

# Check available observations
VO2_data %>%
  count(Condition)

# Fit mixed-effects model
VO2_model <- lmer(
  VO2 ~ Condition + (1 | ID),
  data = VO2_data,
  REML = TRUE
)

# Model results
summary(VO2_model)

# Estimated marginal means and paired contrast
emmeans(VO2_model, pairwise ~ Condition)

# Check model residual normality
shapiro.test(residuals(VO2_model))

# Q-Q plot
qqnorm(residuals(VO2_model))
qqline(residuals(VO2_model), col = "red")


names(Peak)

# Keep valid VO2 measurements
VO2_data <- Peak %>%
  filter(!is.na(VO2), !is.na(Condition))

# Linear mixed-effects model
VO2_model <- lmer(
  VO2 ~ Condition + (1 | ID),
  data = VO2_data
)

# Model results
summary(VO2_model)

# Estimated marginal means
emmeans(VO2_model, ~ Condition)

# Compare CON versus ECC
pairs(emmeans(VO2_model, ~ Condition))


# Shapiro-Wilk test of model residuals
shapiro.test(residuals(VO2_model))

# Q-Q plot
qqnorm(residuals(VO2_model))
qqline(residuals(VO2_model), col = "red")

# Residuals versus fitted values
plot(VO2_model)



# =====================================================
# GXT PEAK - REMAINING PHYSIOLOGICAL VARIABLES
# =====================================================

variables <- c("HR", "RPE", "PO", "Lactate")

for (var in variables) {
  
  cat("\n============================\n")
  cat("VARIABLE:", var, "\n")
  cat("============================\n")
  
  # Select valid observations
  data_var <- Peak[!is.na(Peak[[var]]) &
                     !is.na(Peak$Condition), ]
  
  # Fit linear mixed-effects model
  model <- lmer(
    as.formula(paste(var, "~ Condition + (1 | ID)")),
    data = data_var
  )
  
  # Model summary
  print(summary(model))
  
  # Estimated marginal means
  emm <- emmeans(model, ~ Condition)
  print(emm)
  
  # CON vs ECC comparison
  print(pairs(emm, reverse = TRUE))
  
  # Check residual normality
  print(shapiro.test(residuals(model)))
  
  # Check for singular fit
  print(isSingular(model))
  
  # Q-Q plot of residuals
  qqnorm(residuals(model),
         main = paste("Q-Q Plot:", var))
  qqline(residuals(model), col = "red")
}



# Variables of interest
variables <- c("VO2", "HR", "RPE", "PO", "Lactate")

# Store results
results <- list()

for (var in variables) {
  
  # Keep valid observations
  data_var <- Peak %>%
    filter(!is.na(.data[[var]]),
           !is.na(Condition),
           !is.na(ID))
  
  # Fit LMM
  model <- lmer(
    as.formula(paste(var, "~ Condition + (1 | ID)")),
    data = data_var
  )
  
  # Estimated marginal means
  emm <- emmeans(model, ~ Condition)
  emm_df <- as.data.frame(emm)
  
  # CON - ECC comparison
  comparison <- as.data.frame(
    summary(
      pairs(emm, reverse = TRUE),
      infer = c(TRUE, TRUE)
    )
  )
  
  # Extract results
  results[[var]] <- data.frame(
    Variable = var,
    
    ECC = sprintf(
      "%.2f ± %.2f",
      emm_df$emmean[emm_df$Condition == "ECC"],
      emm_df$SE[emm_df$Condition == "ECC"]
    ),
    
    CON = sprintf(
      "%.2f ± %.2f",
      emm_df$emmean[emm_df$Condition == "CON"],
      emm_df$SE[emm_df$Condition == "CON"]
    ),
    
    Difference = round(comparison$estimate, 2),
    CI_Lower = round(comparison$lower.CL, 2),
    CI_Upper = round(comparison$upper.CL, 2),
    P_value = round(comparison$p.value, 4),
    
    Shapiro_p = round(
      shapiro.test(residuals(model))$p.value, 4
    ),
    
    Singular = isSingular(model)
  )
}

# Combine into one table
GXT_Peak_Results <- bind_rows(results)

# Display
print(GXT_Peak_Results)
View(GXT_Peak_Results)


library(gridExtra)
library(grid)



# Format the table for presentation
Plot_Table <- GXT_Peak_Results %>%
  mutate(
    `CON - ECC` = sprintf("%.2f", Difference),
    `95% CI` = sprintf("[%.2f, %.2f]", CI_Lower, CI_Upper),
    `p-value` = sprintf("%.4f", P_value)
  ) %>%
  select(
    Variable,
    CON,
    ECC,
    `CON - ECC`,
    `95% CI`,
    `p-value`
  )

# Create table
table_plot <- tableGrob(
  Plot_Table,
  rows = NULL,
  theme = ttheme_minimal(
    base_size = 11,
    core = list(
      fg_params = list(hjust = 0.5, x = 0.5)
    ),
    colhead = list(
      fg_params = list(fontface = "bold")
    )
  )
)

# Display in RStudio Plots pane
grid.newpage()

grid.text(
  "GXT Peak Physiological Responses",
  x = 0.5,
  y = 0.92,
  gp = gpar(fontsize = 16, fontface = "bold")
)

grid.draw(
  editGrob(
    table_plot,
    vp = viewport(
      x = 0.5,
      y = 0.48,
      width = 0.95,
      height = 0.72
    )
  )
)



#####Figures#####
# =====================================================
# GXT PEAK FIGURE
# Paired individual data + mean ± SE
# =====================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# Load data
Peak <- read_excel("~/GXT_TTF.xlsx", sheet = "GXT_Peak")

# Set factors
Peak <- Peak %>%
  filter(Condition %in% c("ECC", "CON")) %>%
  mutate(
    ID = factor(ID),
    Condition = factor(Condition, levels = c("ECC", "CON"))
  )

# If needed, set the invalid VO2 reading to NA
# Example only — replace with your actual subject and condition
# Peak <- Peak %>%
#   mutate(
#     VO2 = ifelse(ID == "GMZ_XX" & Condition == "CON", NA, VO2)
#   )

# Variables to plot
plot_data <- Peak %>%
  select(ID, Condition, VO2, HR, RPE, PO, Lactate) %>%
  pivot_longer(
    cols = c(VO2, HR, RPE, PO, Lactate),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value))

# Summary stats for mean ± SE
summary_data <- plot_data %>%
  group_by(Variable, Condition) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    N = sum(!is.na(Value)),
    SE = SD / sqrt(N),
    .groups = "drop"
  )

# Create figure
ggplot(plot_data, aes(x = Condition, y = Value, group = ID)) +
  geom_line(alpha = 0.5) +
  geom_point(size = 2, alpha = 0.8) +
  geom_point(
    data = summary_data,
    aes(x = Condition, y = Mean),
    inherit.aes = FALSE,
    size = 4,
    shape = 18
  ) +
  geom_errorbar(
    data = summary_data,
    aes(
      x = Condition,
      y = Mean,
      ymin = Mean - SE,
      ymax = Mean + SE
    ),
    inherit.aes = FALSE,
    width = 0.08,
    linewidth = 0.8
  ) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 3) +
  theme_classic(base_size = 12) +
  labs(
    title = "Peak Physiological Responses During Eccentric and Concentric GXT",
    x = "Condition",
    y = "Value"
  )



#####BAR PLOTS#####
# Load data
Peak <- read_excel("~/GXT_TTF.xlsx", sheet = "GXT_Peak")

# Set factors
Peak <- Peak %>%
  filter(Condition %in% c("ECC", "CON")) %>%
  mutate(
    ID = factor(ID),
    Condition = factor(Condition, levels = c("ECC", "CON"))
  )

# If needed, set invalid VO2 reading to NA
# Example only:
# Peak <- Peak %>%
#   mutate(
#     VO2 = ifelse(ID == "GMZ_XX" & Condition == "CON", NA, VO2)
#   )

# Convert to long format
plot_data <- Peak %>%
  select(ID, Condition, VO2, HR, RPE, PO, Lactate) %>%
  pivot_longer(
    cols = c(VO2, HR, RPE, PO, Lactate),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value))

# Better labels
plot_data$Variable <- factor(
  plot_data$Variable,
  levels = c("VO2", "HR", "RPE", "PO", "Lactate"),
  labels = c("VO2peak", "HRpeak", "RPEpeak", "Peak Power Output", "Peak Lactate")
)

# Summary stats
summary_data <- plot_data %>%
  group_by(Variable, Condition) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    N = sum(!is.na(Value)),
    SE = SD / sqrt(N),
    .groups = "drop"
  )

# Bar plot
ggplot(summary_data, aes(x = Condition, y = Mean)) +
  geom_col(width = 0.6, alpha = 0.8) +
  geom_errorbar(
    aes(ymin = Mean - SE, ymax = Mean + SE),
    width = 0.15,
    linewidth = 0.8
  ) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 3) +
  theme_classic(base_size = 12) +
  labs(
    title = "Peak Physiological Responses During Eccentric and Concentric GXT",
    x = "Condition",
    y = "Mean ± SE"
  )

#####Try 2 


# =====================================================
# GXT PEAK FIGURE
# Colored bar plots with mean ± SE
# =====================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# Load data
Peak <- read_excel("~/GXT_TTF.xlsx", sheet = "GXT_Peak")

# Set factors
Peak <- Peak %>%
  filter(Condition %in% c("ECC", "CON")) %>%
  mutate(
    ID = factor(ID),
    Condition = factor(Condition, levels = c("ECC", "CON"))
  )

# If needed, set invalid VO2 reading to NA
# Example only:
# Peak <- Peak %>%
#   mutate(
#     VO2 = ifelse(ID == "GMZ_XX" & Condition == "CON", NA, VO2)
#   )

# Convert to long format
plot_data <- Peak %>%
  select(ID, Condition, VO2, HR, RPE, PO, Lactate) %>%
  pivot_longer(
    cols = c(VO2, HR, RPE, PO, Lactate),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value))

# Better labels
plot_data$Variable <- factor(
  plot_data$Variable,
  levels = c("VO2", "HR", "RPE", "PO", "Lactate"),
  labels = c("VO2peak", "HRpeak", "RPEpeak", "Peak Power Output", "Peak Lactate")
)

# Summary stats
summary_data <- plot_data %>%
  group_by(Variable, Condition) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    N = sum(!is.na(Value)),
    SE = SD / sqrt(N),
    .groups = "drop"
  )

# Bar plot
ggplot(summary_data, aes(x = Condition, y = Mean, fill = Condition)) +
  geom_col(width = 0.6, alpha = 0.9) +
  geom_errorbar(
    aes(ymin = Mean - SE, ymax = Mean + SE),
    width = 0.15,
    linewidth = 0.8
  ) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("ECC" = "#4C78A8", "CON" = "#E45756"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Peak Physiological Responses During Eccentric and Concentric GXT",
    x = "Condition",
    y = "Mean ± SE"
  )


#####finallll#####


# =====================================================
# GXT PEAK BAR PLOTS WITH UNITS
# =====================================================

# Define panel labels with units
variable_labels <- c(
  "VO2peak" = "VO2peak (mL/kg/min)",
  "HRpeak" = "HRpeak (bpm)",
  "RPEpeak" = "RPEpeak (Borg 6–20)",
  "Peak Power Output" = "Peak Power Output (W)",
  "Peak Lactate" = "Peak Blood Lactate (mmol/L)"
)

# Create figure
GXT_Peak_Figure <- ggplot(
  summary_data,
  aes(x = Condition, y = Mean, fill = Condition)
) +
  geom_col(
    width = 0.6,
    alpha = 0.9,
    color = "black",
    linewidth = 0.4
  ) +
  geom_errorbar(
    aes(
      ymin = Mean - SE,
      ymax = Mean + SE
    ),
    width = 0.15,
    linewidth = 0.8
  ) +
  facet_wrap(
    ~ Variable,
    scales = "free_y",
    ncol = 3,
    labeller = labeller(Variable = variable_labels)
  ) +
  scale_fill_manual(
    values = c(
      "ECC" = "#4C78A8",
      "CON" = "#E45756"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.25))
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.text = element_text(size = 10)
  ) +
  labs(
    title = "Peak Physiological Responses During GXT",
    x = "Cycling Modality",
    y = "Mean Response"
  )

# Display in Plots pane
print(GXT_Peak_Figure)

# Export high-resolution figure
ggsave(
  "GXT_Peak_Figure.png",
  plot = GXT_Peak_Figure,
  width = 12,
  height = 8,
  dpi = 300
)
