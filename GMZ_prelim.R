library(readxl) ## To load excel sheet
library(dplyr) # Data grammar and manipulation
library(rstatix) # Shapiro Wilk and effect size
library(psych) #descriptives
library(kableExtra) #tables
library(lme4) #linear mixed effects models (LMM)
library(lmerTest) #anova like output for LMM
library(ggplot2) #data visualization
library(ggpubr)#data visualization
library(ggprism)##makes plots look like graphad
library(table1) #for descriptives
library(emmeans)
library(dunn.test)
library(readxl)
library(dplyr)
library(tidyr)
library(rstatix)
library(knitr)
library(kableExtra)


Peak <- read_excel("~/GXT_TTF.xlsx", sheet = "GXT_Peak")
View(Peak)

Peak <- Peak %>%
  mutate(
    ID = factor(ID),
    Condition = factor(
      Condition,
      levels = c("CON", "ECC")
    )
  )

Peak %>%
  count(ID, Condition) %>%
  filter(n > 1)

Peak_wide <- Peak %>%
  pivot_wider(
    names_from = Condition,
    values_from = c(
      PO,
      HR,
      RPE,
      VO2,
      S_BP,
      D_BP,
      Lactate
    ),
    names_glue = "{.value}_{Condition}"
  )

View(Peak_wide)


Peak_wide <- Peak_wide %>%
  mutate(
    PO_diff = PO_ECC - PO_CON,
    HR_diff = HR_ECC - HR_CON,
    RPE_diff = RPE_ECC - RPE_CON,
    VO2_diff = VO2_ECC - VO2_CON,
    S_BP_diff = S_BP_ECC - S_BP_CON,
    D_BP_diff = D_BP_ECC - D_BP_CON,
    Lactate_diff = Lactate_ECC - Lactate_CON
  )

View(Peak_wide)

Peak_wide %>%
  summarise(
    PO_n = sum(complete.cases(PO_CON, PO_ECC)),
    HR_n = sum(complete.cases(HR_CON, HR_ECC)),
    RPE_n = sum(complete.cases(RPE_CON, RPE_ECC)),
    VO2_n = sum(complete.cases(VO2_CON, VO2_ECC)),
    S_BP_n = sum(complete.cases(S_BP_CON, S_BP_ECC)),
    D_BP_n = sum(complete.cases(D_BP_CON, D_BP_ECC)),
    Lactate_n = sum(complete.cases(Lactate_CON, Lactate_ECC))
  )


Peak_normality <- Peak_wide %>%
  select(ID, ends_with("_diff")) %>%
  pivot_longer(
    cols = ends_with("_diff"),
    names_to = "Outcome",
    values_to = "Difference"
  ) %>%
  filter(!is.na(Difference)) %>%
  group_by(Outcome) %>%
  shapiro_test(Difference)

Peak_normality


outcomes <- c(
  "PO",
  "HR",
  "RPE",
  "VO2",
  "S_BP",
  "D_BP",
  "Lactate"
)

run_paired_t <- function(data, outcome) {
  
  con_variable <- paste0(outcome, "_CON")
  ecc_variable <- paste0(outcome, "_ECC")
  
  complete_data <- data %>%
    filter(
      !is.na(.data[[con_variable]]),
      !is.na(.data[[ecc_variable]])
    )
  
  test_result <- t.test(
    complete_data[[ecc_variable]],
    complete_data[[con_variable]],
    paired = TRUE
  )
  
  tibble(
    Outcome = outcome,
    n = nrow(complete_data),
    CON_mean = mean(complete_data[[con_variable]]),
    CON_SD = sd(complete_data[[con_variable]]),
    ECC_mean = mean(complete_data[[ecc_variable]]),
    ECC_SD = sd(complete_data[[ecc_variable]]),
    Mean_difference = mean(
      complete_data[[ecc_variable]] -
        complete_data[[con_variable]]
    ),
    CI_low = test_result$conf.int[1],
    CI_high = test_result$conf.int[2],
    t = unname(test_result$statistic),
    df = unname(test_result$parameter),
    p = test_result$p.value
  )
}

library(purrr)


Peak_t_tests <- map_dfr(
  outcomes,
  ~ run_paired_t(Peak_wide, .x)
)

Peak_t_tests


# ------------------------------------------------------------
# 3. Prepare paired plotting data
# ------------------------------------------------------------

Peak_plot <- Peak %>%
  pivot_longer(
    cols = all_of(outcomes),
    names_to = "Outcome",
    values_to = "Value"
  ) %>%
  group_by(ID, Outcome) %>%
  # Keep only IDs with both CON and ECC for that outcome
  filter(
    sum(!is.na(Value)) == 2,
    n_distinct(Condition[!is.na(Value)]) == 2
  ) %>%
  ungroup()

# Calculate mean and 95% CI using the complete pairs
Peak_summary <- Peak_plot %>%
  group_by(Outcome, Condition) %>%
  summarise(
    n = n(),
    Mean = mean(Value),
    SD = sd(Value),
    SE = SD / sqrt(n),
    CI = qt(
      0.975,
      df = n - 1
    ) * SE,
    .groups = "drop"
  )



# ------------------------------------------------------------
# 4. Create facet labels
# ------------------------------------------------------------

outcome_names <- c(
  PO = "Peak power output (W)",
  HR = "Peak heart rate (bpm)",
  RPE = "Peak RPE (AU)",
  VO2 = "Peak VO2 (mL·kg⁻¹·min⁻¹)",
  S_BP = "Peak systolic BP (mmHg)",
  D_BP = "Peak diastolic BP (mmHg)",
  Lactate = "Peak blood lactate (mmol/L)"
)

Peak_t_tests <- Peak_t_tests %>%
  mutate(
    Outcome_name = unname(outcome_names[Outcome]),
    p_label = ifelse(
      p < 0.001,
      "p < 0.001",
      paste0("p = ", sprintf("%.3f", p))
    ),
    Facet_label = paste0(
      Outcome_name,
      "\nPaired t-test: ",
      p_label,
      "; n = ",
      n
    )
  )

# Add the facet labels to the raw and summary data
Peak_plot <- Peak_plot %>%
  left_join(
    Peak_t_tests %>%
      select(Outcome, Facet_label),
    by = "Outcome"
  )

Peak_summary <- Peak_summary %>%
  left_join(
    Peak_t_tests %>%
      select(Outcome, Facet_label),
    by = "Outcome"
  )




# ------------------------------------------------------------
# 5. Create paired individual-response figure
# ------------------------------------------------------------

GXT_peak_figure <- ggplot(
  Peak_plot,
  aes(
    x = Condition,
    y = Value
  )
) +
  
  # Connect each participant's CON and ECC values
  geom_line(
    aes(group = ID),
    color = "gray70",
    linewidth = 0.7,
    alpha = 0.85
  ) +
  
  # Individual observations
  geom_point(
    aes(color = Condition),
    size = 2.8,
    alpha = 0.9
  ) +
  
  # Group 95% confidence intervals
  geom_errorbar(
    data = Peak_summary,
    aes(
      x = Condition,
      ymin = Mean - CI,
      ymax = Mean + CI
    ),
    inherit.aes = FALSE,
    width = 0.08,
    linewidth = 0.9,
    color = "black"
  ) +
  
  # Group means
  geom_point(
    data = Peak_summary,
    aes(
      x = Condition,
      y = Mean
    ),
    inherit.aes = FALSE,
    shape = 23,
    fill = "white",
    color = "black",
    size = 3.8,
    stroke = 1.2
  ) +
  
  facet_wrap(
    ~ Facet_label,
    scales = "free_y",
    ncol = 4
  ) +
  
  scale_color_manual(
    values = c(
      CON = "#3569A8",
      ECC = "#C4514D"
    )
  ) +
  
  labs(
    title = "Preliminary Peak GXT Responses: CON versus ECC",
    subtitle = paste(
      "Lines connect observations from the same participant.",
      "Diamonds and error bars represent the mean and 95% CI."
    ),
    x = NULL,
    y = NULL,
    caption = paste(
      "P values are unadjusted, two-sided paired t-tests.",
      "Only complete pairs are shown for each outcome."
    )
  ) +
  
  theme_prism(base_size = 12) +
  
  theme(
    legend.position = "none",
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    plot.subtitle = element_text(size = 11),
    plot.caption = element_text(
      size = 9,
      hjust = 0
    ),
    panel.spacing = unit(
      1,
      "lines"
    )
  )

GXT_peak_figure


#####TTF####

TTF_Peak <- read_excel(
  "~/GXT_TTF.xlsx",
  sheet = "TTF_Peak"
)

View(TTF_Peak)

TTF_Peak <- TTF_Peak %>%
  mutate(
    ID = factor(ID),
    Condition = factor(
      Condition,
      levels = c(
        "CON_CON",
        "ECC_CON",
        "ECC_ECC"
      )
    )
  )


TTF_Peak %>%
  count(ID, Condition) %>%
  filter(n > 1)

TTF_Peak %>%
  group_by(Condition) %>%
  summarise(
    Participants = n_distinct(ID),
    PO_n = sum(!is.na(PO)),
    HR_n = sum(!is.na(HR)),
    RPE_n = sum(!is.na(RPE)),
    VO2_n = sum(!is.na(VO2)),
    S_BP_n = sum(!is.na(S_BP)),
    D_BP_n = sum(!is.na(D_BP)),
    Lactate_n = sum(!is.na(Lactate)),
        .groups = "drop"
  )


VO2_model <- lmer(
  VO2 ~ Condition + (1 | ID),
  data = TTF_Peak,
  REML = TRUE,
  na.action = na.exclude
)

summary(VO2_model)

anova

VO2_emmeans <- emmeans(
  VO2_model,
  ~ Condition
)

VO2_pairwise <- pairs(
  VO2_emmeans,
  adjust = "holm"
)

VO2_emmeans
VO2_pairwise

plot(VO2_model)

qqnorm(residuals(VO2_model))
qqline(residuals(VO2_model))


# ============================================================
# TTF PEAK: LMM FIGURE FOR ALL PHYSIOLOGICAL OUTCOMES
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(ggprism)

# ------------------------------------------------------------
# 1. Import and format data
# ------------------------------------------------------------

TTF_Peak <- read_excel(
  "~/GXT_TTF.xlsx",
  sheet = "TTF_Peak"
)

TTF_Peak <- TTF_Peak %>%
  mutate(
    ID = factor(ID),
    Condition = factor(
      Condition,
      levels = c(
        "CON_CON",
        "ECC_CON",
        "ECC_ECC"
      )
    )
  )

# Variables to analyze
outcomes <- c(
  "PO",
  "HR",
  "RPE",
  "VO2",
  "S_BP",
  "D_BP",
  "Lactate"
)

# Readable labels
outcome_names <- c(
  PO = "Peak power output (W)",
  HR = "Peak heart rate (bpm)",
  RPE = "Peak RPE (AU)",
  VO2 = "Peak VO2 (mL·kg⁻¹·min⁻¹)",
  S_BP = "Peak systolic BP (mmHg)",
  D_BP = "Peak diastolic BP (mmHg)",
  Lactate = "Peak blood lactate (mmol/L)"
)


# ------------------------------------------------------------
# 2. Function to fit each LMM
# ------------------------------------------------------------

fit_peak_model <- function(data, outcome) {
  
  model_data <- data %>%
    select(
      ID,
      Condition,
      Value = all_of(outcome)
    ) %>%
    drop_na(Value)
  
  model <- lmer(
    Value ~ Condition + (1 | ID),
    data = model_data,
    REML = TRUE,
    na.action = na.exclude
  )
  
  anova_result <- anova(model)
  
  estimated_means <- emmeans(
    model,
    ~ Condition
  )
  
  list(
    model = model,
    observations = nrow(model_data),
    participants = n_distinct(model_data$ID),
    p = anova_result$`Pr(>F)`[1],
    emmeans = estimated_means
  )
}

# Fit all seven models
TTF_models <- setNames(
  lapply(
    outcomes,
    function(x) fit_peak_model(TTF_Peak, x)
  ),
  outcomes
)


# ------------------------------------------------------------
# 3. Extract omnibus condition tests
# ------------------------------------------------------------

TTF_model_tests <- bind_rows(
  lapply(
    outcomes,
    function(x) {
      
      model_anova <- anova(
        TTF_models[[x]]$model
      )
      
      tibble(
        Outcome = x,
        Observations = TTF_models[[x]]$observations,
        Participants = TTF_models[[x]]$participants,
        F = model_anova$`F value`[1],
        df_num = model_anova$NumDF[1],
        df_den = model_anova$DenDF[1],
        p = model_anova$`Pr(>F)`[1]
      )
    }
  )
)

TTF_model_tests

# ------------------------------------------------------------
# 4. Extract estimated marginal means
# ------------------------------------------------------------

TTF_emmeans <- bind_rows(
  lapply(
    outcomes,
    function(x) {
      
      as.data.frame(
        TTF_models[[x]]$emmeans
      ) %>%
        mutate(
          Outcome = x
        )
    }
  )
)

TTF_emmeans


# ------------------------------------------------------------
# 5. Create panel labels
# ------------------------------------------------------------

TTF_model_tests <- TTF_model_tests %>%
  mutate(
    Outcome_name = unname(
      outcome_names[Outcome]
    ),
    p_label = ifelse(
      p < 0.001,
      "p < 0.001",
      paste0(
        "p = ",
        sprintf("%.3f", p)
      )
    ),
    Facet_label = paste0(
      Outcome_name,
      "\nLMM condition effect: ",
      p_label,
      "; participants = ",
      Participants
    )
  )

TTF_emmeans <- TTF_emmeans %>%
  left_join(
    TTF_model_tests %>%
      select(
        Outcome,
        Facet_label
      ),
    by = "Outcome"
  )

# ------------------------------------------------------------
# 6. Prepare raw data for plotting
# ------------------------------------------------------------

TTF_peak_plot <- TTF_Peak %>%
  pivot_longer(
    cols = all_of(outcomes),
    names_to = "Outcome",
    values_to = "Value"
  ) %>%
  filter(
    !is.na(Value)
  ) %>%
  left_join(
    TTF_model_tests %>%
      select(
        Outcome,
        Facet_label
      ),
    by = "Outcome"
  )


# ------------------------------------------------------------
# 7. Plot individual observations and LMM means
# ------------------------------------------------------------

TTF_peak_figure <- ggplot(
  TTF_peak_plot,
  aes(
    x = Condition,
    y = Value
  )
) +
  
  # Connect available observations from the same participant
  geom_line(
    aes(group = ID),
    color = "gray70",
    linewidth = 0.7,
    alpha = 0.80,
    na.rm = TRUE
  ) +
  
  # Individual observations
  geom_point(
    aes(color = Condition),
    size = 2.7,
    alpha = 0.90
  ) +
  
  # LMM-estimated 95% confidence intervals
  geom_errorbar(
    data = TTF_emmeans,
    aes(
      x = Condition,
      ymin = lower.CL,
      ymax = upper.CL
    ),
    inherit.aes = FALSE,
    width = 0.10,
    linewidth = 0.9,
    color = "black"
  ) +
  
  # LMM-estimated marginal means
  geom_point(
    data = TTF_emmeans,
    aes(
      x = Condition,
      y = emmean
    ),
    inherit.aes = FALSE,
    shape = 23,
    fill = "white",
    color = "black",
    size = 3.8,
    stroke = 1.2
  ) +
  
  facet_wrap(
    ~ Facet_label,
    scales = "free_y",
    ncol = 4
  ) +
  
  scale_x_discrete(
    labels = c(
      CON_CON = "CON–CON",
      ECC_CON = "ECC–CON",
      ECC_ECC = "ECC–ECC"
    )
  ) +
  
  scale_color_manual(
    values = c(
      CON_CON = "#3569A8",
      ECC_CON = "#D18A27",
      ECC_ECC = "#C4514D"
    )
  ) +
  
  labs(
    title = "Preliminary Peak TTF Responses Across Conditions",
    subtitle = paste(
      "Lines connect available observations from the same participant.",
      "Diamonds and error bars represent LMM-estimated means and 95% CIs."
    ),
    x = NULL,
    y = NULL,
    caption = paste(
      "P values represent the omnibus condition effect from",
      "random-intercept linear mixed-effects models."
    )
  ) +
  
  theme_prism(base_size = 12) +
  
  theme(
    legend.position = "none",
    strip.text = element_text(
      face = "bold",
      size = 9.5
    ),
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    ),
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    plot.subtitle = element_text(
      size = 11
    ),
    plot.caption = element_text(
      size = 9,
      hjust = 0
    )
  )

TTF_peak_figure