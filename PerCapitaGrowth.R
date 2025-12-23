##Code for the 90 % confidence envelope of breaks==========================================
#And stack graphs of mean KDE over time, 90 % confidence
#envelope of per capita growth over time, and mean per capita growth vs.
#mean KDE density (plot C). Plot C completes steps 1 and 2
#in the main manuscript to estimate logistic growth functions and 
#transitions to functions with higher K. All plots were annotated and
#placed into the directory FinalSIFIgs.

##Load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(cowplot)

#load data
growth <- read.csv("data/Percapita2/AllCasesFoodPro200sims2.csv")

# Ensure PeriodID exists
if ("PeriodID" %in% names(growth)) {
  growth$PeriodID <- as.factor(growth$PeriodID)
} else {
  stop("Variable 'PeriodID' does not exist in the dataset.")
}

# Identify simulation columns
pop_cols <- grep("^V[0-9]+$", names(growth), value = TRUE)

# --------------------------------------------------------------
# 2) Compute per-simulation growth matrix for each region
# --------------------------------------------------------------
compute_region_gmatrix <- function(df_region, pop_cols, arrange_ascending = FALSE) {
  if(arrange_ascending) df_region <- df_region %>% arrange(calBP)
  else df_region <- df_region %>% arrange(desc(calBP))
  
  mat <- as.matrix(df_region[, pop_cols, drop = FALSE])
  ntime <- nrow(mat)
  nsim  <- ncol(mat)
  if (ntime < 2) stop("Region has fewer than 2 time points.")
  
  g_mat <- matrix(NA_real_, nrow = ntime - 1, ncol = nsim)
  
  for (j in seq_len(nsim)) {
    vec <- as.numeric(mat[, j])
    g_j <- log(vec[-1] / vec[-ntime])
    g_j[is.infinite(g_j)] <- NA_real_
    g_mat[, j] <- g_j
  }
  
  list(calBP_g = df_region$calBP[-ntime], g_mat = g_mat)
}

# --------------------------------------------------------------
# 3) Summaries: mean, 10%, 90% growth
# --------------------------------------------------------------
summarize_gmatrix <- function(calBP_g, g_mat) {
  
  out <- tibble(
    calBP = calBP_g,
    mean_g = rowMeans(g_mat, na.rm = TRUE),
    lo10 = apply(g_mat, 1, quantile, probs = 0.10, na.rm = TRUE),
    hi90 = apply(g_mat, 1, quantile, probs = 0.90, na.rm = TRUE)
  )
  
  # *** APPLY 0.7 CAP HERE ***
  out <- out %>%
    mutate(
      mean_g = pmin(mean_g, 0.6),
      lo10   = pmin(lo10, 0.6),
      hi90   = pmin(hi90, 0.6)
    )
  
  out <- out %>%
    mutate(
      mean_g = pmax(mean_g, -0.5),
      lo10   = pmax(lo10, -0.5),
      hi90   = pmax(hi90, -0.5)
    )
  
  return(out)
}

# --------------------------------------------------------------
# 4) Growth summary for all regions
# --------------------------------------------------------------
regions <- unique(growth$region_id2)

growth_summary <- map_dfr(regions, function(rid) {
  df_region <- filter(growth, region_id2 == rid)
  gm <- compute_region_gmatrix(df_region, pop_cols, arrange_ascending = FALSE)
  summarize_gmatrix(gm$calBP_g, gm$g_mat) %>%
    mutate(region_id2 = rid) %>%
    select(region_id2, everything())
})

# --------------------------------------------------------------
# 5) KDE summary (mean + standardized)
# --------------------------------------------------------------
kde_summary <- map_dfr(regions, function(rid) {
  df_region <- filter(growth, region_id2 == rid) %>% arrange(desc(calBP))
  
  mat <- as.matrix(df_region[, pop_cols, drop = FALSE])
  
  mean_kde <- rowMeans(mat, na.rm = TRUE)
  lo10_kde <- apply(mat, 1, quantile, probs = 0.10, na.rm = TRUE)
  hi90_kde <- apply(mat, 1, quantile, probs = 0.90, na.rm = TRUE)
  
  mn <- min(mean_kde, na.rm = TRUE)
  mx <- max(mean_kde, na.rm = TRUE)
  rng <- mx - mn
  
  if (rng == 0 || is.na(rng)) {
    StKDE <- StKDE_lo <- StKDE_hi <- rep(0, length(mean_kde))
  } else {
    StKDE <- (mean_kde - mn)/rng
    StKDE_lo <- (lo10_kde - mn)/rng
    StKDE_hi <- (hi90_kde - mn)/rng
  }
  
  tibble(
    region_id2 = rid,
    calBP = df_region$calBP,
    mean_kde, lo10_kde, hi90_kde,
    StKDE, StKDE_lo, StKDE_hi,
    PeriodID = df_region$PeriodID
  )
})

write.table(kde_summary, file = "data/Percapita2/KDESummarytest.csv", sep = ",", col.names=NA)

# --------------------------------------------------------------
# 6) Merge growth and KDE summaries
# --------------------------------------------------------------
combined_df <- growth_summary %>%
  left_join(kde_summary, by = c("region_id2", "calBP"))

# --------------------------------------------------------------
# 7) Output directory
# --------------------------------------------------------------
out_dir <- "TestingSIFigs"
if (!dir.exists(out_dir)) dir.create(out_dir)

# --------------------------------------------------------------
# 8) Theme
# --------------------------------------------------------------
my_theme <- theme_bw() +
  theme(
    axis.text.x  = element_text(size=28, colour="black"),
    axis.text.y  = element_text(size=28, colour="black"),
    axis.title.x = element_text(size=24),
    axis.title.y = element_text(size=24),
    plot.title   = element_text(size=22, face="bold"),
    legend.title = element_text(size=20),
    legend.text  = element_text(size=18)
  )

# --------------------------------------------------------------
# STORAGE FOR EQUILIBRIUM SUMMARY FOR ALL REGIONS
# --------------------------------------------------------------
all_equilibria <- list()

# --------------------------------------------------------------
# 9) Produce per-region plots
# --------------------------------------------------------------
plot_regions <- unique(combined_df$region_id2)

for (reg in plot_regions) {
  
  df <- filter(combined_df, region_id2 == reg)
  if (nrow(df) < 2) next
  
  df <- df %>% arrange(desc(calBP))
  
  # ------------------------------------------------------------
  # (A) Standardized KDE ONLY
  # ------------------------------------------------------------
  pA <- ggplot(df, aes(x = calBP)) +
    geom_ribbon(aes(ymin = StKDE_lo, ymax = StKDE_hi),
                fill="gray70", alpha=0.6) +
    geom_line(aes(y = StKDE), color="black", linewidth=1.1) +
    geom_point(aes(y = StKDE, color = PeriodID), size = 2) +
    scale_x_reverse() +
    labs(title=paste("(A) Standardized KDE Density VS. Time –", reg),
         x="Years calBP", y="Standardized KDE density", color="Period") +
    my_theme
  
  # ------------------------------------------------------------
  # (B) Per-capita growth vs time
  # ------------------------------------------------------------
  pB <- ggplot(df, aes(x = calBP)) +
    geom_ribbon(aes(ymin=lo10, ymax=hi90), fill="gray70", alpha=0.6) +
    geom_line(aes(y=mean_g), linewidth=0.9) +
    geom_point(aes(y=mean_g, color=PeriodID), size=2) +
    scale_x_reverse() +
    geom_hline(yintercept=0) +
    labs(title=paste("(B) KDE Per capita Growth VS. Time –", reg),
         x="Years calBP", y="KDE per capita growth", color="Period") +
    my_theme
  
  # ------------------------------------------------------------
  # (C) Growth vs Standardized KDE — WITH FILTERED EQUILIBRIA
  # ------------------------------------------------------------
  df2 <- df %>% select(StKDE, mean_g)
  
  # Find exact zeros
  zero_exact <- df2 %>% 
    mutate(idx = row_number()) %>% 
    filter(mean_g == 0)
  
  # Interpolated zero crossings
  interp <- list()
  for (i in seq_len(nrow(df2)-1)) {
    y1 <- df2$mean_g[i];  y2 <- df2$mean_g[i+1]
    x1 <- df2$StKDE[i];   x2 <- df2$StKDE[i+1]
    
    if (y1 == 0)
      interp[[length(interp)+1]] <- data.frame(StKDE=x1, mean_g=0, idx=i)
    
    if (y2 == 0)
      interp[[length(interp)+1]] <- data.frame(StKDE=x2, mean_g=0, idx=i+1)
    
    if (sign(y1) == sign(y2)) next
    if ((y2-y1)==0) next
    
    x0 <- x1 - y1*(x2 - x1)/(y2 - y1)
    interp[[length(interp)+1]] <- data.frame(StKDE=x0, mean_g=0, idx=i+0.5)
  }
  
  zero_interp <- bind_rows(interp)
  zero_points_raw <- bind_rows(zero_exact, zero_interp) %>% arrange(idx)
  
  # ---- FILTER: require 2 consecutive positive g prior ----
  keep_equilibria <- function(df2, zeros) {
    if (nrow(zeros) == 0) return(zeros)
    out <- list()
    for (i in seq_len(nrow(zeros))) {
      idx <- zeros$idx[i]
      prev1 <- floor(idx) - 1
      prev2 <- floor(idx) - 2
      if (prev2 >= 1 &&
          df2$mean_g[prev1] > 0 &&
          df2$mean_g[prev2] > 0) {
        out[[length(out)+1]] <- zeros[i,]
      }
    }
    bind_rows(out)
  }
  
  zero_points <- keep_equilibria(df2, zero_points_raw)
  
  # anchor lines
  anchor <- data.frame(StKDE=0, mean_g=0.3)
  
  if (nrow(zero_points) > 0) {
    segments <- bind_rows(lapply(seq_len(nrow(zero_points)), function(i) {
      data.frame(
        StKDE = c(anchor$StKDE, zero_points$StKDE[i]),
        mean_g = c(anchor$mean_g, zero_points$mean_g[i]),
        seg_id = i
      )
    }))
  } else {
    segments <- tibble(StKDE=numeric(), mean_g=numeric())
  }
  
  # Save equilibria
  if (nrow(zero_points) > 0) {
    all_equilibria[[reg]] <- zero_points %>% mutate(region = reg)
  }
  
  pC <- ggplot(df, aes(x = StKDE, y = mean_g)) +
    geom_path(linewidth=1.0) +
    geom_point(aes(color=PeriodID), size=2) +
    geom_hline(yintercept=0) +
    {if (nrow(segments)>0)
      geom_line(data=segments, aes(x=StKDE,y=mean_g,group=seg_id),
                color="darkgray", linetype="dashed", linewidth=1)} +
    {if (nrow(zero_points)>0)
      geom_point(data=zero_points, aes(x=StKDE,y=mean_g),
                 color="darkgray", size=4)} +
    labs(title=paste("(C) KDE Per Capita Growth vs. Density –", reg),
         x="Standardized KDE density", y="KDE per capita growth", color="Period") +
    my_theme
  
  # ------------------------------------------------------------
  # Stack A, B, C
  # ------------------------------------------------------------
  final_plot <- plot_grid(pA, pB, pC, ncol=1, align="v", axis="l")
  
  ggsave(file.path(out_dir, paste0(reg, ".pdf")),
         final_plot, width=11.55, height=14, units="in")
}



