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


####GXT_DURING####

# ============================================================
# GXT DURING: ALL AVAILABLE STAGE RESPONSES
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggprism)

# ------------------------------------------------------------
# 1. Import data
# ------------------------------------------------------------

GXT_During <- read_excel(
  "~/GXT_TTF.xlsx",
  sheet = "GXT_During"
)

stage_levels <- c(
  "Baseline",
  "One",
  "Two",
  "Three",
  "Four",
  "Five",
  "Six",
  "Seven",
  "Eight"
)

outcomes <- c(
  "PO",
  "HR",
  "RPE",
  "VO2",
  "S_BP",
  "D_BP",
  "Lactate"
)

outcome_labels <- c(
  "Peak power output (W)",
  "Heart rate (bpm)",
  "RPE (AU)",
  "VO2 (mL·kg⁻¹·min⁻¹)",
  "Systolic BP (mmHg)",
  "Diastolic BP (mmHg)",
  "Blood lactate (mmol/L)"
)

# ------------------------------------------------------------
# 2. Format conditions and stages
# ------------------------------------------------------------

GXT_During <- GXT_During %>%
  filter(
    !is.na(ID),
    !is.na(Condition),
    !is.na(Stage)
  ) %>%
  mutate(
    ID = factor(ID),
    
    Condition = factor(
      Condition,
      levels = c("CON", "ECC")
    ),
    
    Stage = factor(
      trimws(as.character(Stage)),
      levels = stage_levels,
      ordered = TRUE
    ),
    
    Stage_number = as.numeric(Stage) - 1
  )

GXT_During %>%
  filter(
    S_BP > 300
  ) %>%
  select(
    ID,
    Condition,
    Stage,
    S_BP
  )


# ------------------------------------------------------------
# 3. Convert to long format
# ------------------------------------------------------------

GXT_long <- GXT_During %>%
  pivot_longer(
    cols = all_of(outcomes),
    names_to = "Outcome",
    values_to = "Value"
  ) %>%
  filter(
    !is.na(Value)
  ) %>%
  mutate(
    Outcome = factor(
      Outcome,
      levels = outcomes,
      labels = outcome_labels
    )
  )


# ------------------------------------------------------------
# 4. Stage-level descriptive statistics
# ------------------------------------------------------------

GXT_summary <- GXT_long %>%
  group_by(
    Outcome,
    Condition,
    Stage_number
  ) %>%
  summarise(
    n = n_distinct(ID),
    Mean = mean(Value),
    SD = sd(Value),
    SE = SD / sqrt(n),
    CI = ifelse(
      n > 1,
      qt(0.975, df = n - 1) * SE,
      NA_real_
    ),
    .groups = "drop"
  )

GXT_summary


GXT_stage_counts <- GXT_long %>%
  group_by(
    Outcome,
    Condition,
    Stage_number
  ) %>%
  summarise(
    n = n_distinct(ID),
    .groups = "drop"
  ) %>%
  arrange(
    Outcome,
    Condition,
    Stage_number
  )

GXT_stage_counts


# ------------------------------------------------------------
# 5. Plot individual and group trajectories
# ------------------------------------------------------------

GXT_during_figure <- ggplot(
  GXT_long,
  aes(
    x = Stage_number,
    y = Value
  )
) +
  
  # Individual participant trajectories
  geom_line(
    aes(
      group = interaction(ID, Condition),
      color = Condition
    ),
    linewidth = 0.45,
    alpha = 0.20
  ) +
  
  geom_point(
    aes(color = Condition),
    size = 1.3,
    alpha = 0.25
  ) +
  
  # Group mean trajectories
  geom_line(
    data = GXT_summary,
    aes(
      x = Stage_number,
      y = Mean,
      color = Condition,
      linetype = Condition,
      group = Condition
    ),
    inherit.aes = FALSE,
    linewidth = 1.25
  ) +
  
  # Group means
  geom_point(
    data = GXT_summary,
    aes(
      x = Stage_number,
      y = Mean,
      color = Condition,
      shape = Condition
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  
  # Group 95% confidence intervals
  geom_errorbar(
    data = GXT_summary,
    aes(
      x = Stage_number,
      ymin = Mean - CI,
      ymax = Mean + CI,
      color = Condition
    ),
    inherit.aes = FALSE,
    width = 0.10,
    linewidth = 0.65,
    na.rm = TRUE
  ) +
  
  facet_wrap(
    ~ Outcome,
    scales = "free_y",
    ncol = 4
  ) +
  
  scale_x_continuous(
    breaks = 0:8,
    labels = stage_levels
  ) +
  
  scale_color_manual(
    values = c(
      CON = "#3569A8",
      ECC = "#C4514D"
    )
  ) +
  
  scale_shape_manual(
    values = c(
      CON = 16,
      ECC = 17
    )
  ) +
  
  scale_linetype_manual(
    values = c(
      CON = "solid",
      ECC = "dashed"
    )
  ) +
  
  labs(
    title = "Preliminary Physiological Responses During the GXT",
    subtitle = paste(
      "Thin lines represent individual participants;",
      "bold lines represent observed group means."
    ),
    x = "GXT stage",
    y = NULL,
    color = "Condition",
    shape = "Condition",
    linetype = "Condition",
    caption = paste(
      "All available observations are shown.",
      "Later stages contain fewer participants because GXT duration differed."
    )
  ) +
  
  theme_prism(base_size = 11) +
  
  theme(
    legend.position = "bottom",
    strip.text = element_text(
      face = "bold",
      size = 9
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    plot.caption = element_text(
      size = 9,
      hjust = 0
    )
  )

GXT_during_figure



#####GXT during analysis####

library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(knitr)
library(kableExtra)

# ------------------------------------------------------------
# 1. Define stages and outcomes
# ------------------------------------------------------------

analysis_stages <- c(
  "Baseline",
  "One",
  "Two",
  "Three",
  "Four",
  "Five",
  "Six"
)

outcomes <- c(
  "PO",
  "HR",
  "RPE",
  "VO2",
  "S_BP",
  "D_BP",
  "Lactate"
)

GXT_analysis <- GXT_During %>%
  filter(
    !is.na(ID),
    !is.na(Condition),
    !is.na(Stage),
    as.character(Stage) %in% analysis_stages
  ) %>%
  mutate(
    ID = factor(ID),
    Condition = factor(
      Condition,
      levels = c("CON", "ECC")
    ),
    Stage = factor(
      as.character(Stage),
      levels = analysis_stages,
      ordered = FALSE
    )
  ) %>%
  droplevels()


# ------------------------------------------------------------
# 2. Fit Condition × Stage LMM
# ------------------------------------------------------------

run_stage_lmm <- function(data, outcome) {
  
  model_data <- data %>%
    select(
      ID,
      Condition,
      Stage,
      Value = all_of(outcome)
    ) %>%
    filter(
      !is.na(Value)
    ) %>%
    droplevels()
  
  model <- lmer(
    Value ~ Condition * Stage + (1 | ID),
    data = model_data,
    REML = TRUE,
    na.action = na.exclude
  )
  
  # Overall model effects
  omnibus <- as.data.frame(
    anova(
      model,
      type = 3,
      ddf = "Satterthwaite"
    )
  ) %>%
    tibble::rownames_to_column(
      var = "Effect"
    )
  
  # Estimated means for CON and ECC within every stage
  stage_emmeans <- emmeans(
    model,
    ~ Condition | Stage
  )
  
  # ECC minus CON at every stage
  stage_contrasts <- contrast(
    stage_emmeans,
    method = list(
      "ECC - CON" = c(-1, 1)
    ),
    adjust = "none"
  ) %>%
    as.data.frame()
  
  # Number of participants contributing to each condition-stage
  sample_sizes <- model_data %>%
    group_by(
      Stage,
      Condition
    ) %>%
    summarise(
      n = n_distinct(ID),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = Condition,
      values_from = n,
      names_prefix = "n_",
      values_fill = 0
    )
  
  stage_contrasts <- stage_contrasts %>%
    left_join(
      sample_sizes,
      by = "Stage"
    )
  
  list(
    model = model,
    omnibus = omnibus,
    emmeans = stage_emmeans,
    contrasts = stage_contrasts
  )
}

# ------------------------------------------------------------
# 3. Run every outcome
# ------------------------------------------------------------

GXT_stage_models <- setNames(
  lapply(
    outcomes,
    function(x) {
      run_stage_lmm(
        GXT_analysis,
        x
      )
    }
  ),
  outcomes
)


# ------------------------------------------------------------
# 4. Extract omnibus tests
# ------------------------------------------------------------

GXT_omnibus_results <- bind_rows(
  lapply(
    outcomes,
    function(x) {
      
      GXT_stage_models[[x]]$omnibus %>%
        mutate(
          Outcome = x,
          .before = 1
        )
    }
  )
)

GXT_omnibus_results %>%
  select(
    Outcome,
    Effect,
    NumDF,
    DenDF,
    `F value`,
    `Pr(>F)`
  )

# ------------------------------------------------------------
# 5. Extract stage-specific ECC-CON comparisons
# ------------------------------------------------------------

GXT_stage_comparisons <- bind_rows(
  lapply(
    outcomes,
    function(x) {
      
      GXT_stage_models[[x]]$contrasts %>%
        mutate(
          Outcome = x,
          .before = 1
        )
    }
  )
)

GXT_stage_comparisons <- GXT_stage_comparisons %>%
  group_by(Outcome) %>%
  mutate(
    p_Holm_across_stages = p.adjust(
      p.value,
      method = "holm"
    )
  ) %>%
  ungroup() %>%
  mutate(
    across(
      c(
        estimate,
        SE,
        df,
        t.ratio
      ),
      ~ round(.x, 2)
    ),
    p_raw = round(
      p.value,
      3
    ),
    p_Holm = round(
      p_Holm_across_stages,
      3
    )
  )

GXT_stage_comparisons %>%
  select(
    Outcome,
    Stage,
    n_CON,
    n_ECC,
    contrast,
    estimate,
    SE,
    df,
    t.ratio,
    p_raw,
    p_Holm
  )

GXT_stage_comparisons %>%
  filter(
    Outcome == "VO2"
  )

GXT_stage_comparisons %>%
  filter(
    Outcome == "Lactate"
  )

GXT_stage_comparisons %>%
  filter(
    Outcome == "HR"
  )

GXT_stage_comparisons %>%
  filter(
    Outcome %in% c(
      "PO",
      "HR",
      "RPE",
      "VO2",
      "Lactate"
    )
  ) %>%
  select(
    Outcome,
    Stage,
    n_CON,
    n_ECC,
    estimate,
    SE,
    df,
    t.ratio,
    p_raw,
    p_Holm
  ) %>%
  arrange(
    factor(
      Outcome,
      levels = c(
        "PO",
        "HR",
        "RPE",
        "VO2",
        "S_BP",
        "D_BP",
        "Lactate"
      )
    ),
    Stage
  ) %>%
  kbl(
    caption = "Model-estimated ECC versus CON comparisons during the GXT"
  ) %>%
  kable_classic(
    full_width = FALSE,
    html_font = "Cambria"
  )
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


library(ggpubr)

# ============================================================
# PAIRWISE COMPARISONS FROM THE LMMs
# ============================================================

TTF_pairwise <- bind_rows(
  lapply(
    outcomes,
    function(x) {
      
      pairwise_result <- pairs(
        TTF_models[[x]]$emmeans,
        adjust = "holm"
      )
      
      as.data.frame(pairwise_result) %>%
        mutate(
          Outcome = x
        )
    }
  )
)

TTF_pairwise


TTF_pairwise <- TTF_pairwise %>%
  separate(
    contrast,
    into = c("group1", "group2"),
    sep = " - ",
    remove = FALSE
  ) %>%
  mutate(
    p_label = case_when(
      p.value < 0.001 ~ "Holm p < 0.001",
      TRUE ~ paste0(
        "Holm p = ",
        sprintf("%.3f", p.value)
      )
    )
  )

TTF_pairwise <- TTF_pairwise %>%
  left_join(
    TTF_model_tests %>%
      select(
        Outcome,
        Facet_label
      ),
    by = "Outcome"
  )


TTF_plot_ranges <- TTF_peak_plot %>%
  group_by(
    Outcome,
    Facet_label
  ) %>%
  summarise(
    y_min = min(Value, na.rm = TRUE),
    y_max = max(Value, na.rm = TRUE),
    y_range = y_max - y_min,
    .groups = "drop"
  ) %>%
  mutate(
    y_range = ifelse(
      y_range == 0,
      abs(y_max),
      y_range
    )
  )

TTF_plot_ranges <- TTF_Peak %>%
  pivot_longer(
    cols = all_of(outcomes),
    names_to = "Outcome",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value)) %>%
  group_by(Outcome) %>%
  summarise(
    y_min = min(Value),
    y_max = max(Value),
    y_range = y_max - y_min,
    .groups = "drop"
  ) %>%
  mutate(
    y_range = ifelse(
      y_range == 0,
      abs(y_max),
      y_range
    )
  )

TTF_plot_ranges

TTF_pairwise <- TTF_pairwise %>%
  left_join(
    TTF_plot_ranges,
    by = "Outcome"
  ) %>%
  group_by(Outcome) %>%
  arrange(Outcome, contrast) %>%
  mutate(
    y.position = y_max +
      c(0.10, 0.23, 0.36) * y_range
  ) %>%
  ungroup()

TTF_pairwise

TTF_peak_pairwise_figure <- TTF_peak_figure +
  stat_pvalue_manual(
    TTF_pairwise,
    label = "p_label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.4,
    size = 3,
    hide.ns = FALSE
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0.08, 0.42)
    )
  )

TTF_peak_pairwise_figure


####TTF DUring ####

library(readxl)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(knitr)
library(kableExtra)

# ------------------------------------------------------------
# 1. Import and format
# ------------------------------------------------------------

TTF_During <- read_excel(
  "~/GXT_TTF.xlsx",
  sheet = "TTF_During"
)

all_stage_levels <- c(
  "Baseline",
  "One",
  "Two",
  "Three",
  "Four",
  "Five",
  "Six",
  "Seven",
  "Eight",
  "Nine"
)

TTF_During <- TTF_During %>%
  filter(
    !is.na(ID),
    !is.na(Condition),
    !is.na(Stage)
  ) %>%
  mutate(
    ID = factor(ID),
    Condition = factor(
      Condition,
      levels = c(
        "CON_CON",
        "ECC_CON",
        "ECC_ECC"
      )
    ),
    Stage = factor(
      trimws(as.character(Stage)),
      levels = all_stage_levels,
      ordered = FALSE
    )
  )

# ------------------------------------------------------------
# 2. Common stages for three-condition inference
# ------------------------------------------------------------

TTF_stage_counts <- TTF_During %>%
  pivot_longer(
    cols = c(
      HR,
      RPE,
      VO2,
      S_BP,
      D_BP,
      Lactate
    ),
    names_to = "Outcome",
    values_to = "Value"
  ) %>%
  filter(
    !is.na(Value)
  ) %>%
  count(
    Outcome,
    Stage,
    Condition,
    name = "n"
  ) %>%
  pivot_wider(
    names_from = Condition,
    values_from = n,
    values_fill = 0
  ) %>%
  arrange(
    Outcome,
    Stage
  )

TTF_stage_counts

print(
  TTF_stage_counts,
  n = Inf
)

Outcome ~ Condition * Stage + (1 | ID)

common_stages <- c(
  "Baseline",
  "One",
  "Two"
)

# ------------------------------------------------------------
# 2. Retain stages represented in all three conditions
# ------------------------------------------------------------

common_stages <- c(
  "Baseline",
  "One",
  "Two"
)

TTF_common <- TTF_During %>%
  filter(
    as.character(Stage) %in% common_stages
  ) %>%
  mutate(
    Stage = factor(
      as.character(Stage),
      levels = common_stages,
      ordered = FALSE
    )
  ) %>%
  droplevels()

analysis_outcomes <- c(
  "HR",
  "RPE",
  "VO2",
  "S_BP",
  "D_BP",
  "Lactate"
)


# ------------------------------------------------------------
# 3. LMM function
# ------------------------------------------------------------

run_TTF_stage_model <- function(data, outcome) {
  
  model_data <- data %>%
    select(
      ID,
      Condition,
      Stage,
      Value = all_of(outcome)
    ) %>%
    filter(!is.na(Value)) %>%
    droplevels()
  
  model <- lmer(
    Value ~ Condition * Stage + (1 | ID),
    data = model_data,
    REML = TRUE,
    na.action = na.exclude,
    control = lmerControl(
      optimizer = "bobyqa"
    )
  )
  
  # Overall condition, stage, and interaction tests
  omnibus <- anova(
    model,
    type = 3,
    ddf = "Satterthwaite"
  ) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(
      var = "Effect"
    )
  
  # Model-estimated condition means within each stage
  estimated_means <- emmeans(
    model,
    ~ Condition | Stage
  )
  
  # Three condition comparisons within each stage
  comparisons <- pairs(
    estimated_means,
    adjust = "none"
  ) %>%
    as.data.frame()
  
  # Actual sample size for each condition-stage
  sample_sizes <- model_data %>%
    count(
      Stage,
      Condition,
      name = "n"
    ) %>%
    pivot_wider(
      names_from = Condition,
      values_from = n,
      names_prefix = "n_",
      values_fill = 0
    )
  
  comparisons <- comparisons %>%
    left_join(
      sample_sizes,
      by = "Stage"
    )
  
  list(
    model = model,
    omnibus = omnibus,
    emmeans = estimated_means,
    comparisons = comparisons
  )
}

# ------------------------------------------------------------
# 4. Run all outcome models
# ------------------------------------------------------------

TTF_stage_models <- setNames(
  lapply(
    analysis_outcomes,
    function(x) {
      run_TTF_stage_model(
        TTF_common,
        x
      )
    }
  ),
  analysis_outcomes
)

names(TTF_stage_models)


# ------------------------------------------------------------
# 5. Extract overall LMM tests
# ------------------------------------------------------------

TTF_omnibus_results <- bind_rows(
  lapply(
    analysis_outcomes,
    function(x) {
      
      TTF_stage_models[[x]]$omnibus %>%
        mutate(
          Outcome = x,
          .before = 1
        )
    }
  )
)

TTF_omnibus_results_clean <- TTF_omnibus_results %>%
  transmute(
    Outcome,
    Effect,
    df_num = NumDF,
    df_den = DenDF,
    F = `F value`,
    p = `Pr(>F)`
  ) %>%
  mutate(
    across(
      c(df_den, F),
      ~ round(.x, 2)
    ),
    p = round(p, 3)
  )

TTF_omnibus_results_clean


TTF_omnibus_results_clean %>%
  kbl(
    caption = paste(
      "Linear mixed-effects models of physiological",
      "responses during the common TTF stages"
    )
  ) %>%
  kable_classic(
    full_width = FALSE,
    html_font = "Cambria"
  )


# ------------------------------------------------------------
# 6. Extract model-based pairwise comparisons
# ------------------------------------------------------------

TTF_stage_comparisons <- bind_rows(
  lapply(
    analysis_outcomes,
    function(x) {
      
      TTF_stage_models[[x]]$comparisons %>%
        mutate(
          Outcome = x,
          .before = 1
        )
    }
  )
)


# Holm adjustment across the three condition comparisons
# separately within each stage
TTF_stage_comparisons <- TTF_stage_comparisons %>%
  group_by(
    Outcome,
    Stage
  ) %>%
  mutate(
    p_Holm_within_stage = p.adjust(
      p.value,
      method = "holm"
    )
  ) %>%
  ungroup()

# More conservative Holm adjustment across all nine
# condition-stage comparisons within each outcome
TTF_stage_comparisons <- TTF_stage_comparisons %>%
  group_by(Outcome) %>%
  mutate(
    p_Holm_across_stages = p.adjust(
      p.value,
      method = "holm"
    )
  ) %>%
  ungroup()


TTF_stage_comparisons_clean <- TTF_stage_comparisons %>%
  transmute(
    Outcome,
    Stage,
    n_CON_CON,
    n_ECC_CON,
    n_ECC_ECC,
    Comparison = contrast,
    Estimate = estimate,
    SE,
    df,
    t = t.ratio,
    p_raw = p.value,
    p_Holm_within_stage,
    p_Holm_across_stages
  ) %>%
  mutate(
    across(
      c(
        Estimate,
        SE,
        df,
        t
      ),
      ~ round(.x, 2)
    ),
    across(
      c(
        p_raw,
        p_Holm_within_stage,
        p_Holm_across_stages
      ),
      ~ round(.x, 3)
    )
  ) %>%
  arrange(
    factor(
      Outcome,
      levels = analysis_outcomes
    ),
    Stage,
    Comparison
  )

print(
  TTF_stage_comparisons_clean,
  n = Inf
)

TTF_stage_comparisons_clean %>%
  kbl(
    caption = paste(
      "Model-estimated condition comparisons",
      "within the common TTF stages"
    )
  ) %>%
  kable_classic(
    full_width = FALSE,
    html_font = "Cambria"
  )







# ============================================================
# DESCRIPTIVE TTF DURING FIGURE
# No inferential analysis
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggprism)

# ------------------------------------------------------------
# 1. Import data
# ------------------------------------------------------------

TTF_During <- read_excel(
  "~/GXT_TTF.xlsx",
  sheet = "TTF_During"
)

stage_levels <- c(
  "Baseline",
  "One",
  "Two",
  "Three",
  "Four",
  "Five",
  "Six",
  "Seven",
  "Eight",
  "Nine"
)

outcomes <- c(
  "PO",
  "HR",
  "RPE",
  "VO2",
  "S_BP",
  "D_BP",
  "Lactate"
)

outcome_labels <- c(
  "Power output (W)",
  "Heart rate (bpm)",
  "RPE (AU)",
  "VO2 (mL·kg⁻¹·min⁻¹)",
  "Systolic BP (mmHg)",
  "Diastolic BP (mmHg)",
  "Blood lactate (mmol/L)"
)


# ------------------------------------------------------------
# 2. Format variables
# ------------------------------------------------------------

TTF_During <- TTF_During %>%
  filter(
    !is.na(ID),
    !is.na(Condition),
    !is.na(Stage)
  ) %>%
  mutate(
    ID = factor(ID),
    
    Condition = factor(
      Condition,
      levels = c(
        "CON_CON",
        "ECC_CON",
        "ECC_ECC"
      )
    ),
    
    Stage = factor(
      trimws(as.character(Stage)),
      levels = stage_levels,
      ordered = TRUE
    ),
    
    Stage_number = as.numeric(Stage) - 1
  )

# ------------------------------------------------------------
# 3. Create long plotting dataset
# ------------------------------------------------------------

TTF_plot_data <- TTF_During %>%
  pivot_longer(
    cols = all_of(outcomes),
    names_to = "Outcome",
    values_to = "Value"
  ) %>%
  filter(
    !is.na(Value)
  ) %>%
  mutate(
    Outcome = factor(
      Outcome,
      levels = outcomes,
      labels = outcome_labels
    )
  )


# ------------------------------------------------------------
# 4. Calculate observed means
# ------------------------------------------------------------

TTF_plot_summary <- TTF_plot_data %>%
  group_by(
    Outcome,
    Condition,
    Stage_number
  ) %>%
  summarise(
    n = n_distinct(ID),
    Mean = mean(
      Value,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

TTF_plot_summary


# ------------------------------------------------------------
# 5. Plot every available observation
# ------------------------------------------------------------

TTF_during_descriptive_figure <- ggplot(
  TTF_plot_data,
  aes(
    x = Stage_number,
    y = Value
  )
) +
  
  # Individual participant trajectories
  geom_line(
    aes(
      group = interaction(ID, Condition),
      color = Condition
    ),
    linewidth = 0.55,
    alpha = 0.25
  ) +
  
  # Individual participant observations
  geom_point(
    aes(color = Condition),
    size = 1.5,
    alpha = 0.35
  ) +
  
  # Observed group mean trajectories
  geom_line(
    data = TTF_plot_summary,
    aes(
      x = Stage_number,
      y = Mean,
      color = Condition,
      linetype = Condition,
      group = Condition
    ),
    inherit.aes = FALSE,
    linewidth = 1.4
  ) +
  
  # Observed group means
  geom_point(
    data = TTF_plot_summary,
    aes(
      x = Stage_number,
      y = Mean,
      color = Condition,
      shape = Condition
    ),
    inherit.aes = FALSE,
    size = 3.2,
    stroke = 1
  ) +
  
  facet_wrap(
    ~ Outcome,
    scales = "free_y",
    ncol = 4
  ) +
  
  scale_x_continuous(
    breaks = 0:9,
    labels = stage_levels
  ) +
  
  scale_color_manual(
    values = c(
      CON_CON = "#3569A8",
      ECC_CON = "#D18A27",
      ECC_ECC = "#C4514D"
    ),
    labels = c(
      CON_CON = "CON–CON",
      ECC_CON = "ECC–CON",
      ECC_ECC = "ECC–ECC"
    )
  ) +
  
  scale_shape_manual(
    values = c(
      CON_CON = 16,
      ECC_CON = 17,
      ECC_ECC = 15
    ),
    labels = c(
      CON_CON = "CON–CON",
      ECC_CON = "ECC–CON",
      ECC_ECC = "ECC–ECC"
    )
  ) +
  
  scale_linetype_manual(
    values = c(
      CON_CON = "solid",
      ECC_CON = "dashed",
      ECC_ECC = "dotdash"
    ),
    labels = c(
      CON_CON = "CON–CON",
      ECC_CON = "ECC–CON",
      ECC_ECC = "ECC–ECC"
    )
  ) +
  
  labs(
    title = "Physiological Responses Throughout the TTF Trials",
    subtitle = paste(
      "Thin lines represent individual participants;",
      "thick lines and large points represent observed means."
    ),
    x = "TTF stage",
    y = NULL,
    color = "Condition",
    shape = "Condition",
    linetype = "Condition",
    caption = paste(
      "All available observations are displayed.",
      "Trajectories end when the participant reached task failure;",
      "later-stage means therefore contain fewer participants."
    )
  ) +
  
  theme_prism(base_size = 11) +
  
  theme(
    legend.position = "bottom",
    
    strip.text = element_text(
      face = "bold",
      size = 9.5
    ),
    
    axis.text.x = element_text(
      angle = 45,
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

TTF_during_descriptive_figure