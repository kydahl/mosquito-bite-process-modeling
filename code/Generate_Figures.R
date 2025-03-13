# Code to generate the figures in ``Once bitten, twice shy: A modeling framework for incorporating heterogeneous mosquito biting into transmission models''

# Load libraries ----
library(tidyverse)
library(latex2exp)
library(cols4all)
library(cowplot)
library(expm)
library(matlib)
library(scales)
library(ggh4x)

# Load data ----
Full_df = read_rds("data/GCD_R0_data.rds") %>% 
  mutate(varied_parameter = if_else(
    is.na(varied_parameter),
    "none",
    varied_parameter)) %>% 
  # Set up labels
  mutate(Type = if_else(
    `Model type` == "Mechanistic",
    paste0(`Model type`, " (", varied_parameter,")"),
    `Model type`))

# Dictionary for labeling parameters nicely
param_table = tibble(
  Symbol = c("$p_Q$", "$\\lambda_Q$", "$p_L$", "$\\lambda_L$", "$p_P$", "$\\lambda_P$", "$p_G$", "$\\lambda_G$", "$\\sigma$", "$b$"),
  Description = c(
    "Probability of progressing from seeking to landing",
    "Exit rate from seeking stage (per minute)", 
    "Probability of progressing from landing to probing",
    "Exit rate from landing stage (per minute)", 
    "Probability of progressing from probing to ingesting",
    "Exit rate from probing stage (per minute)",
    "Probability of progressing from ingesting to ovipositing",
    "Exit rate from ingestion stage (per minute)",
    "Probability of seeking a new host given feeding failure",
    "Biting rate (exponential model)"
    
  ),
  Label = c("Seeking success", "Seeking rate", "Landing success", "Landing rate", "Probing success", "Probing rate", "Ingesting success", "Ingesting rate", "Persistence probability", "Biting rate (exponential)"), 
  Type = c("B Probabilities", "A Rates", "B Probabilities", "A Rates", "B Probabilities", "A Rates", "B Probabilities", "A Rates", "B Probabilities", "A Rates"),
  Prefix = c("Seeking", "Seeking", "Landing", "Landing", "Probing", "Probing", "Ingesting",  "Ingesting", "Persistence", "Exponential"), 
  short_label = c("pQ", "lQ", "pL", "lL", "pP", "lP", "pG", "lG", "sigma", "b")
)

# Set varied_parameter levels for mechanistic model
Full_df$varied_parameter = factor(Full_df$varied_parameter, 
                                  levels = c("none", "lQ", "pP", "pG"))

# Labels and ordering for the model types
type_order = c("Standard", "Exponential", "Empirical", "Phenomenological", "Mechanistic")
Full_df$`Model type` = factor(Full_df$`Model type`, 
                              levels = type_order)

# Define functions ----

# Function: probability distribution function for a phase-type distribution
PH_pdf <- function(x, A_matrix, v_alpha) {
  A_dim = dim(A_matrix)[1]
  if (is.null(A_dim)) {
    v_ones = 1
    pdf_val = exp(A_matrix * x) * (-A_matrix)
  } else {
    v_ones = matrix(rep(1, A_dim), ncol = 1)
    
    pdf_val = t(v_alpha) %*% expm(A_matrix * x) %*% (-A_matrix %*% v_ones)
  }
  return(pdf_val)
}

# Function: mean of a phase-type distribution
PH_mean <- function(A_matrix, v_alpha) {
  A_dim = dim(A_matrix)[1]
  if (is.null(A_dim)) {
    v_ones = 1
    mean_val = (-1 / A_matrix[[1]])
  } else {
    v_ones = matrix(rep(1, A_dim), ncol = 1)
    
    mean_val = t(v_ones) %*% inv(-t(A_matrix)) %*% v_alpha
  }
  return(mean_val)
}

# Helper function: place legends in empty facets of plot grids
# Code by: Artem Sokolov, found here: https://stackoverflow.com/questions/54438495/shift-legend-into-empty-facets-of-a-faceted-plot-in-ggplot2
shift_legend <- function(p) {
  pnls <- cowplot::plot_to_gtable(p) %>%
    gtable::gtable_filter("panel") %>%
    with(setNames(grobs, layout$name)) %>%
    purrr::keep(~ identical(.x, zeroGrob()))
  
  if (length(pnls) == 0) stop("No empty facets in the plot")
  
  lemon::reposition_legend(p, "center", panel = names(pnls))
}

# Plot parameters ----

# Colors for the model types
type_colors = c("black", c4a("brewer.dark2",4))

# Figure 4. Pdfs of distributions ----

# Plot the pdf of each model type for a couple 
# values of theta: 1/4 days, 1/2 days, 2 days
theta_vals = c((1/4) * 1440, (1/2) * 1440, 1440, 2 * 1440)

# Resolution and maximum value to use for drawing the pdf
pdf_max = 8 * 1440
resolution = 10000

# For mechanistic model, we can't guarantee we'll hit the theta values exactly, so just choose the closest ones
Mech_df_filtered <- Full_df %>% 
  filter(`Model type` %in% c("Mechanistic")) %>% 
  filter(varied_parameter == "lQ") %>% # the curves for lQ, pP, pG look almost identical, so just look at one
  group_by(`Model type`, varied_parameter) %>% 
  filter(theta %in% sapply(theta_vals, function(x) theta[which.min(abs(theta - x))])) %>%
  ungroup()

# Assemble data for plotting
Figure4_df <- Full_df %>% 
  # Remove the Mechanistic rows, which are handled separately
  filter(!(`Model type` %in% c("Mechanistic"))) %>% 
  # Combine Standard and Exponential since they are equivalent in this case
  mutate(Type = if_else(
    Type %in% c("Standard", "Exponential"),
    "Standard / Exponential",
    Type)
  ) %>% 
  # Just keep chosen theta values
  filter(theta %in% theta_vals) %>% 
  # Add the filtered Mechanistic rows
  bind_rows(Mech_df_filtered) %>% 
  # Assign the closest theta to each row (to relate Mechanistic with the rest)
  group_by(Type) %>% 
  mutate(closest_theta = theta_vals[sapply(theta, function(x) which.min(abs(x - theta_vals)))]) %>% 
  ungroup() %>% 
  # Select only relevant rows
  dplyr::select(Type, theta, closest_theta, v_alpha, A_matrix) %>% 
  mutate(theta_label = case_when(
    closest_theta == 360 ~ "A",
    closest_theta == 720 ~ "B",
    closest_theta == 1440 ~ "C",
    closest_theta == 2880 ~ "D"
  )) %>%
  # Add in values to plug into the pdf
  cross_join(tibble(x = seq(0, pdf_max, length.out = resolution)))  %>% 
  filter(case_when(
    closest_theta == 360 ~ x < 0.75 * 1440,
    closest_theta == 720 ~ x < 1.5 * 1440,
    closest_theta == 1440 ~ x < 3 * 1440,
    closest_theta == 2880 ~ x < 6 * 1440
  )) %>%
  # Calculate values of the pdf
  rowwise() %>% 
  mutate(pdf_val = PH_pdf(x, A_matrix, v_alpha))

# Define levels for model types
Figure4_df$Type = factor(
  Figure4_df$Type,
  levels = c("Empirical", "Phenomenological", "Mechanistic (lQ)","Standard / Exponential"))
# Make nice labels for model types
Figure4_labels = c(expression("Standard / Exponential"), expression("Empirical"), expression("Phenomenological"), 
                   expression("Mechanistic " (lambda[Q]))
)
# Define colors for model types
Fig4_color_vals = c("Standard / Exponential" = "black",
                    "Empirical" = c4a("brewer.dark2", 3)[1],
                    "Phenomenological" = c4a("brewer.dark2", 3)[2],
                    "Mechanistic (lQ)" = c4a("brewer.dark2", 3)[3]#
)

# set up labels for the subplots
facet_labels <- Figure4_df %>%
  group_by(closest_theta) %>% 
  filter(x > 0) %>% 
  mutate(label = LETTERS[1:n()],  # Assign "A", "B", "C", ...
         x_pos = 360*min(x, na.rm = TRUE)/1440,  # Align left
         y_pos = max(pdf_val, na.rm = TRUE) * .95) %>%  # Slightly above max y
  distinct(closest_theta, x_pos, y_pos, theta_label)

# Figure 1: pdfs of model types for some values of theta
Figure4 = Figure4_df %>%
  filter(closest_theta < 2 * 1440) %>% 
  ggplot(aes(x = 24 * x / (1440), y = pdf_val, color = Type)) +
  # Plot grey dotted line showing the mean
  geom_vline(
    aes(xintercept = 24 * closest_theta/1440),
    lwd = 0.5, color = "grey", lty = 2
  ) +
  # Plot pdfs
  geom_line(
    lwd = 0.5,
    alpha = 0.75
  ) + 
  # Label subplots A, B, C, ...
  geom_text(data = facet_labels %>% filter(closest_theta < 2 * 1440),
            aes(x = 0.01, y = y_pos, label = theta_label, group = theta_label),
            color = "black", size = 4,
            hjust ="inward",
  ) +
  # Wrap by theta value
  facet_wrap( ~ closest_theta,
              nrow = 1,
              scales = "free") +
  scale_x_continuous(
    name = "Gonotrophic cycle duration [Hours]",
    limits =  c(0, NA),
    expand = c(0.0,0),
    breaks = seq(0, 72, by = 6)
  ) +
  scale_y_continuous(
    name = "Density",
    limits = c(0, NA),
    expand = c(0.01,0),
  ) +
  scale_color_manual(
    name = "Model type:",
    values = Fig4_color_vals,
    breaks = unique(Figure4_df$Type),
    labels = Figure4_labels
  ) +
  guides(
    # Put the color legend at the top of the plot
    color = guide_legend(
      position = "top",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  theme_half_open(10) +
  theme(
    strip.background = element_rect(color = "white", fill = "white"),
    strip.text = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

Figure4
ggsave("figures/Figure4.pdf", Figure4, width = 6.5, height = 1.75 * 9/6.5, units = "in")

# Figure 5. R0 vs. standard biting rates ----

# Plot R0 as a function of the standard biting rate (1/theta) for each model type
# Distinguish types by color and linetype
# Make sure the labels are nice

# Colors for model types
Fig5_color_vals = c("Standard" = "black",
                    "Exponential" = "black",
                    "Empirical" = c4a("brewer.dark2", 3)[1],
                    "Phenomenological" = c4a("brewer.dark2", 3)[2],
                    "Mechanistic~(lambda[Q])" = c4a("brewer.dark2", 3)[3],
                    "Mechanistic~(p[P])" = c4a("brewer.dark2", 3)[3],
                    "Mechanistic~(p[G])" = c4a("brewer.dark2", 3)[3])
# Linetypes for model types
Fig5_lty_vals = c("Standard" = 1,
                  "Exponential" = 2,
                  "Empirical" = 1,
                  "Phenomenological" = 1,
                  "Mechanistic~(lambda[Q])" = 2,
                  "Mechanistic~(p[P])" = 3,
                  "Mechanistic~(p[G])" = 4)

# Assemble data for figure 2
Figure5_df = Full_df %>%
  # Make labels nice for mechanistic model
  mutate(
    Type = case_when(
      Type == "Mechanistic (lQ)" ~ "Mechanistic~(lambda[Q])",
      Type == "Mechanistic (pP)" ~ "Mechanistic~(p[P])",
      Type == "Mechanistic (pG)" ~ "Mechanistic~(p[G])",
      TRUE ~ Type
    ),
    sbr = 1440/theta  # standard biting rate
  )

# Collect values where R0 first exceeds or less deceeds one to label on the x-axis
Figure5_ticks <- Figure5_df %>% 
  filter(sbr < 2.01) %>% # only keep values up to 2 bites per day
  # Get x-coordinate where R0 first exceeds/is less than one
  group_by(Type) %>% 
  arrange(sbr) %>% 
  summarise(
    first_R0_greater_1 = sbr[R0 > 1][1],
    first_R0_less_1 = max(sbr[R0 > 1])
  ) %>% 
  unique()

# Nice labels for model types
Figure5_labels = c(expression("Standard"), expression("Exponential"), expression("Empirical"), expression("Phenomenological"), 
                   expression("Mechanistic " (lambda[Q])), expression("Mechanistic " (p[P])),expression("Mechanistic " (p[G]))
)

# Arrows showing the direction of increase of sbr and R0 as a mechanistic parameter is increased
Figure5_arrows = Figure5_df %>%
  filter(`Model type` == "Mechanistic") %>% 
  pivot_longer(cols = lQ:sigma) %>% 
  filter(varied_parameter == name) %>% 
  group_by(varied_parameter) %>%
  arrange(value) %>% 
  mutate(sbr = 1440 / theta) %>% 
  filter(between(sbr, 1.36, 1.42)) %>% 
  group_by(varied_parameter) %>% 
  mutate(
    x_first = min(sbr),
    x_last = max(sbr)
  ) %>% 
  group_by(varied_parameter) %>% 
  filter(sbr %in% c(x_first, x_last)) %>% 
  mutate(
    R0_first = R0[sbr == min(sbr)],
    R0_last = R0[sbr == max(sbr)]
  ) %>% ungroup() %>% 
  select(Type, x_first, x_last, R0_first, R0_last) %>% distinct()

# Figure 2: R0 vs GCD for all model types
Figure5 <- Figure5_df %>% 
  filter(sbr < 2.01) %>% # only keep values up to 2 bites per day
  ggplot(aes(color = Type, lty = Type)) +
  # Grey line for R0 = 1
  geom_hline(aes(yintercept = 1), color = "grey", lwd = 1) +
  # R0-GCD curves
  geom_line(aes(x = 1440 / theta, y = R0),
            lwd = 0.75) +
  # Add ticks below the x-axis for R0 = 1 crossings
  geom_rug(
    data = Figure5_ticks,
    aes(x = first_R0_greater_1),
    sides = "b", size = 0.75, outside = TRUE,
    length = unit(0.3, "in"),
    show.legend = F
  ) +
  # Arrows showing direction of increasing parameter values
  geom_segment(
    data = Figure5_arrows,
    aes(x = x_first, y = R0_first, xend = x_last, yend = R0_last, color = Type),
    lwd = 0, 
    # alpha = 0,
    arrow = arrow(length = unit(0.2, "inches"), ends = "last", type = "closed"),
    show.legend = F,
    inherit.aes = F
  ) +
  scale_x_continuous(
    name = TeX("Standard biting rate [Days$^{-1}$]"),
    expand = c(0, 0, 0, 0.01)
  ) +
  scale_y_continuous(
    name = TeX("Basic reproduction number \\, [$R_0$]"),
    limits = c(0, NA),
    expand = c(0,0)
  ) +
  scale_linetype_manual(
    name = "Model type:",
    values = Fig5_lty_vals,
    breaks = unique(Figure5_df$Type),
    labels = Figure5_labels
  ) +
  scale_color_manual(
    name = "Model type:",
    values = Fig5_color_vals,
    breaks = unique(Figure5_df$Type),
    labels = Figure5_labels
  ) +
  coord_cartesian(
    xlim = c(0, 1.8),
    clip = "off" # needed to let ticks get plotted out of the axes
  ) + 
  theme_half_open(11) + 
  theme(
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.key.width = unit(0.4, "in")
  ) +
  guides(
    color = guide_legend(
      position = "top",
      direction = "horizontal",
      nrow = 2,
      byrow = T
    ),
    linetype = guide_legend(
      position = "top",
      direction = "horizontal",
      nrow = 2,
      byrow = T
    )
  )

# Save Figure 5
ggsave("figures/Figure5.pdf", Figure5, width = 7.5, height = 3.25 * 9/6.5, units = "in")


# Table 4: Characteristics of R0 curves ----
R0_characteristics_table <- Full_df %>%
  filter(between(theta, 1, 20*1440)) %>% # remove unrealistically short GCD (less than 1 second)
  mutate(sbr = 1440/theta) %>%  # standard biting rate
  arrange(sbr) %>% 
  mutate(
    Type = case_when(
      Type == "Mechanistic (lQ)" ~ "Mechanistic~(lambda[Q])",
      Type == "Mechanistic (pP)" ~ "Mechanistic~(p[P])",
      Type == "Mechanistic (pG)" ~ "Mechanistic~(p[G])",
      TRUE ~ Type
    )
  ) %>% 
  # Get x-coordinate where R0 first exceeds/is less than one
  group_by(Type) %>% 
  summarise(
    max_sbr = max(sbr),
    crit_min_sbr = sbr[R0 > 1][1],
    max_test = length(sbr[R0>1]),
    temp_max = max(sbr[R0 > 1], na.rm = T),
    crit_max_sbr = ifelse(temp_max > 0.99*max_sbr | max_test == 0, NA, temp_max),
    max_R0 = max(R0)
  ) %>% 
  select(-c(max_sbr, max_test, temp_max)) %>% 
  unique()

write_csv(R0_characteristics_table, "data/Table4.csv")

# Figure 6. R0 vs. mechanistic parameters ----
# Load in all data for mechanistic model
Mech_df <- read_rds("data/Mechanistic_results.rds")

# Distinguish parameters by color -- different from 2.
nice_mech_labels = data.frame(
  pQ = "Seeking~success*','~p[Q]",
  pL = "Landing~success*','~p[L]",
  pP = "Probing~success*','~p[P]",
  pG = "Ingesting~success*','~p[G]",
  sigma = "Persistence~probability*','~sigma",
  lQ = "Seeking~rate*','~lambda[Q]",
  lL = "Landing~rate*','~lambda[L]",
  lP = "Probing~rate*','~lambda[P]",
  lG = "Ingesting~rate*','~lambda[G]"
) %>% pivot_longer(everything(), values_to = "nice_labels")

# Assemble data for figure 6
Figure6_df <- Mech_df %>% 
  filter(parameter_type == "varied") %>% 
  select(-parameter_type) %>% 
  pivot_longer(lQ:sigma) %>% 
  filter(name == varied_parameter) %>% 
  mutate(parameter_type = case_when(
    name %in% c("lQ", "lL", "lP", "lG") ~ "rate",
    name %in% c("pQ", "pL", "pP", "pG", "sigma") ~ "probability"
  )) %>% 
  # Reduce the range for rates
  filter(!(parameter_type == "rate" & value >= (1/60))) %>% 
  right_join(nice_mech_labels)

# Assign levels for parameters names
Figure6_df$name <- factor(
  Figure6_df$name,
  levels = c("pQ", "pL", "pP", "pG", "sigma", "lQ", "lL", "lP", "lG"))

# Assign nice labels
Figure6_df$nice_labels <- factor(
  Figure6_df$nice_labels,
  levels = nice_mech_labels$nice_labels)

# Figure 6
Figure6 <- Figure6_df %>% 
  ggplot(aes(x = value, y = R0, linetype = mosquito_type, color = parameter_type)) +
  geom_hline(aes(yintercept = 1), color = "grey") +
  geom_line(lwd = 0.5) +
  facet_wrap(
    ~ nice_labels,
    ncol = 5,
    labeller = labeller(nice_labels = label_parsed),
    scales = "free_x"
  ) +
  scale_x_continuous(
    name = "",
    breaks = waiver(),
    n.breaks = 5,
    labels = function(x) sub("\\.?0+$", "", format(x, nsmall = 2)),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = TeX("Basic reproduction number \\, [$R_0$]"),
    labels = function(x) sub("\\.?0+$", "", format(x, nsmall = 2)),
    expand = c(0,0)
  ) +
  scale_color_manual(
    name = "Parameter type:",
    values = c4a("parks.saguaro",2)
  ) +
  scale_linetype_discrete(
    name = "Mosquito type:",
    labels = c("Flighty", "Persistent")
  ) +
  theme_half_open(font_size = 8) +
  guides(
    color = guide_none()
  ) +
  theme(
    strip.background = element_rect(color = "white", fill = "white"),
    legend.key.width = unit(0.4, "in")
  )

shift_legend(Figure6)
ggsave("figures/Figure6.pdf", shift_legend(Figure6), width = 6.5, height = 2.25 * 9/6.5, units = "in")

# Figure 7. PRCCs of R0 against mechanistic parameters ----

# Set this to TRUE to calculate the PRCC values from the LHS samples produced by the Julia code
calculate_PRCCs_bool = TRUE

if (calculate_PRCCs_bool) {
  # Load in data
  LHS_data = read_csv("data/LHS_samples.csv.gz") %>%
    filter(type == "max")

  rank_data = LHS_data %>%
    group_by(type) %>%
    mutate(across(lQ:R0, ~ rank(.x)))

  PRCC_data <- rank_data %>%
    pivot_longer(cols = lQ:dummy, names_to = "input", values_to = "input_value") %>%
    pivot_longer(cols = GCD:R0, names_to = "output", values_to = "output_value") %>%
    ungroup() %>%
    # group_by(type, input, output) %>%
    summarise(
      PRCC = cor(input_value, output_value),
      .by = c(type, input, output)
    )

  # Save the final PRCC results
  write_csv(PRCC_data, "data/PRCC_data.csv")
} else {
  # Load in final PRCC results
  PRCC_data <- read_csv("data/PRCC_data.csv")
}

# Assemble data for Figure 7
plot_data <- PRCC_data %>% 
  group_by(type, output) %>% 
  # Add in nice labels
  left_join(
    param_table %>% 
      rename(input = short_label)
  ) %>%
  left_join(rename(nice_mech_labels, input = name), by = "input")  %>% 
  mutate(
    output_label = case_when(
      output == "R0" ~ "Basic reproduction number",
      output == "GCD" ~ "Gonotrophic cycle duration",
      output == "N_offspring" ~ "Basic offspring number",
    ),
    type_label = case_when(
      type == "flighty" ~ "Flighty",
      type == "persistent" ~ "Persistent",
      type == "max" ~ "Maximum variation",
    ),
    dummy_min = -abs(PRCC[input == "dummy"]),
    dummy_max = abs(PRCC[input == "dummy"])
  ) %>% 
  filter(!is.na(output), input != "dummy") %>% 
  mutate(input_num = as.numeric(factor(input)))
  
# Assign levels for parameter set types
plot_data$type = factor(plot_data$type, levels = c("flighty", "persistent", "max"))
# Assign levels for parameters
plot_data$input = factor(plot_data$input, levels = c(
  "pQ", "pL", "pP", "pG", "sigma", "lQ", "lL", "lP", "lG", "dummy"
))
# Assign levels for outputs
plot_data$output_label = factor(plot_data$output_label, levels = c(
  "Gonotrophic cycle duration", "Basic offspring number","Basic reproduction number"
))
# Assign levels for type labels
plot_data$type_label = factor(plot_data$type_label, levels = c(
  "Flighty","Persistent","Maximum variation"
))
# Assign labels for parameters
plot_data$Label = factor(plot_data$Label, levels = rev(c(
  c("Seeking success",
    "Landing success",  "Probing success", "Ingesting success"),
  "Persistence probability", 
  "Seeking rate", "Landing rate", "Probing rate", "Ingesting rate", "Dummy variable"
)))
# Order the labels
plot_data$nice_labels <- factor(
  plot_data$nice_labels,
  levels = c(rev(nice_mech_labels$nice_labels), "Dummy~variable"))

Figure7 <- plot_data %>% 
  # Just keep maximum variation for manuscript
  filter(type == "max") %>%
  arrange(input) %>% 
  filter(output %in% c("N_offspring","R0")) %>%
  group_by(type, output) %>% 
  mutate(star_flag = abs(PRCC) < abs(dummy_min),
         star_xpos = PRCC + sign(PRCC) * 0.025) %>% 
  # Plot
  ggplot() +
  geom_col(aes(y = nice_labels, x = PRCC, fill = output_label), position = "dodge", width = 0.75) +
  # Add light grey lines to divide up categories
  geom_hline(yintercept = seq(1.5, length(levels(plot_data$input)) - 0.5, by = 1), 
             color = "grey60", linetype = "dashed", linewidth = 0.125) +
  # Add zero line
  geom_vline(xintercept = 0, color = "black", lwd = 0.25) +
  scale_fill_manual(
    name = "",
    values = c(c4a("met.juarez",3))
  ) +
  # Add "n.s." for flagged 'insignificant' bars
  geom_text(
    data = . %>% filter(star_flag),
    aes(y = nice_labels, x = star_xpos, label = "n.s."),
    position = position_dodge(width = 0.75),
    size = 2, vjust = -0.55
  ) +
  scale_y_discrete(
    name = "",
    labels = function(x) parse(text = x)
  ) +
  scale_x_continuous(
    TeX("Partial Rank Correlation Coefficient"),
    breaks = seq(-1,1,by = 0.25),
    expand = c(0.05,0)
  ) +
  theme_half_open(11) +
  theme(
    strip.background = element_rect(color = "white", fill = "white"),
    legend.key.width = unit(0.25, "in"),
    legend.key.height = unit(0.03125, "in"),
    legend.position = "top",
    legend.direction = "horizontal"
  )

Figure7

ggsave("figures/Figure7.pdf", Figure7, width = 6.5, height = 2.25 * 9/6.5, units = "in")