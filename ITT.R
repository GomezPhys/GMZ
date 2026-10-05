library(dplyr)
library(readxl) ## To load excel sheet

MVC <- read.csv(
  "~/GMZ10_ECC_ECC_POST_ITT.csv",
  check.names = FALSE
)

View(MVC)

# Convert required columns to numeric
MVC_data <- MVC %>%
  mutate(
    `End Pnt 0` = as.numeric(`End Pnt 0`),
    `Time(Seconds)` = as.numeric(`Time(Seconds)`),
    `Torque(Foot-Pounds)` = as.numeric(`Torque(Foot-Pounds)`)
  )

# Find the row containing the highest torque in each repetition
MVC_peaks <- MVC_data %>%
  filter(
    `End Pnt 0` %in% c(0, 2, 4),
    !is.na(`Torque(Foot-Pounds)`)
  ) %>%
  group_by(`End Pnt 0`) %>%
  slice_max(
    order_by = `Torque(Foot-Pounds)`,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(
    Repetition = case_when(
      `End Pnt 0` == 0 ~ "Rep 1",
      `End Pnt 0` == 2 ~ "Rep 2",
      `End Pnt 0` == 4 ~ "Rep 3"
    )
  ) %>%
  select(
    Repetition,
    `End Pnt 0`,
    Peak_Torque_ft_lb = `Torque(Foot-Pounds)`,
    Time_of_Peak_sec = `Time(Seconds)`
  )

# Formatted copy for viewing
MVC_peaks_display <- MVC_peaks %>%
  mutate(
    Peak_Torque_ft_lb = sprintf("%.4f", Peak_Torque_ft_lb),
    Time_of_Peak_sec = sprintf("%.3f", Time_of_Peak_sec)
  )

print(MVC_peaks_display)

# Find the highest MVC across all three repetitions
Peak_MVC <- MVC_peaks %>%
  slice_max(
    order_by = Peak_Torque_ft_lb,
    n = 1,
    with_ties = FALSE
  )

Peak_MVC_display <- Peak_MVC %>%
  mutate(
    Peak_Torque_ft_lb = sprintf("%.4f", Peak_Torque_ft_lb),
    Time_of_Peak_sec = sprintf("%.3f", Time_of_Peak_sec)
  )

cat("\nOverall peak MVC:\n")
print(Peak_MVC_display)


# ============================================================
# PLOT TORQUE DEVELOPMENT FOR EACH ITT REPETITION
# ============================================================

library(ggplot2)

# Prepare the three repetition traces
ITT_plot_data <- MVC_data %>%
  filter(`End Pnt 0` %in% c(0, 2, 4)) %>%
  mutate(
    Repetition = case_when(
      `End Pnt 0` == 0 ~ "Rep 1",
      `End Pnt 0` == 2 ~ "Rep 2",
      `End Pnt 0` == 4 ~ "Rep 3"
    ),
    Repetition = factor(
      Repetition,
      levels = c("Rep 1", "Rep 2", "Rep 3")
    )
  ) %>%
  arrange(Repetition, `Time(Seconds)`)

# Create the figure
ITT_torque_plot <- ggplot(
  ITT_plot_data,
  aes(
    x = `Time(Seconds)`,
    y = `Torque(Foot-Pounds)`
  )
) +
  geom_line(
    color = "#1F5A94",
    linewidth = 0.6
  ) +
  geom_hline(
    yintercept = 0,
    color = "gray50",
    linetype = "dashed",
    linewidth = 0.4
  ) +
  facet_wrap(
    ~ Repetition,
    ncol = 1,
    scales = "free_x"
  ) +
  labs(
    title = "Torque Development During ITT Trials",
    x = "Time (seconds)",
    y = "Torque (ft-lb)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(
      fill = "gray90",
      color = "black"
    ),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

# Display the figure
print(ITT_torque_plot)