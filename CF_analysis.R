# CF_analysis.R
# Purpose: Causal forest analysis of HS employment on medium-run outcomes (ELS:2002 public-use)
# Outcomes: Y_bachelors_plus (BA or higher by 2012), Y_income (annual earnings 2011)
# Sections:
#   1. BY treatment  - any work in 10th grade (D_by), sample: cleaner_by.dta
#   2. F1 treatment  - any work in 12th grade (D_f1), sample: cleaner_f1.dta
#   3. MACF          - work intensity in 12th grade (W_f1), sample: cleaner_f1.dta
# Output: output/by_cf/ and output/f1_cf/

# ---- USER: SET THIS PATH BEFORE RUNNING ----
wd <- "C:\\Users\\pardonic\\Documents\\HS employment"
# --------------------------------------------

set.seed(1)
setwd(wd)

dir.create("output/by_cf",  showWarnings = FALSE, recursive = TRUE)
dir.create("output/f1_cf",  showWarnings = FALSE, recursive = TRUE)

sink("temp/CF_analysis_log.txt", split = TRUE)
msg_con <- file("temp/CF_analysis_errors.txt", open = "wt")
sink(msg_con, type = "message")
cat("CF_analysis.R started:", format(Sys.time()), "\n\n")

# ==============================================================================
# 0. Setup — packages, data, feature matrices, helper functions
# ==============================================================================

pkgs <- c("haven", "grf", "dplyr", "tidyr", "ggplot2", "stringr", "weights",
          "RColorBrewer", "survey", "sandwich", "lmtest", "patchwork")
pkgs_missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(pkgs_missing) > 0) {
  install.packages(pkgs_missing, repos = "https://cloud.r-project.org")
}
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- Load data ----
df_by <- read_dta("temp/cleaner_by.dta")
df_f1 <- read_dta("temp/cleaner_f1.dta")

# Extract Stata variable labels
var_labels <- sapply(df_by, function(x) {
  lbl <- attr(x, "label")
  if (length(lbl) == 0 || all(lbl == "")) NULL else lbl[1]
})
var_labels <- unlist(var_labels[!sapply(var_labels, is.null)])

df_by <- haven::zap_labels(df_by) #save base year data frame with labels fixed
df_f1 <- haven::zap_labels(df_f1) #save 12th grade data frame w labels fixed

get_label <- function(v) {
  lbl <- var_labels[v]
  unname(ifelse(is.na(lbl), v, lbl))
}

cat("BY sample loaded:", nrow(df_by),
    "| D_by treated =", sum(df_by$D_by),
    "/ control =", sum(df_by$D_by == 0), "\n")
cat("F1 sample loaded:", nrow(df_f1),
    "| D_f1 treated =", sum(df_f1$D_f1),
    "/ control =", sum(df_f1$D_f1 == 0), "\n\n")
#display in log sample size and number who worked

# ---- Feature matrix definition ----
vl_cat_prefixes <- c("byrace", "byfcomp", "bysctrl", "byregurb", "byocc30") #categorical variables which have this prefix
dummy_pattern   <- paste0("^(", paste(vl_cat_prefixes, collapse = "|"), ")_\\d+$") #get the pattern wiht _##
dummy_cols      <- grep(dummy_pattern, names(df_by), value = TRUE) 
dummy_cols      <- c(dummy_cols, "bymale", "bystlang") #the dummy variables that do not follow that pattern, bound to those that do 

vl_cont <- c("byses2", "bytxcstd", "byscsaf2",
             "bystprep", "bytstrel", "byxtracu", "f1age",
             "bypared", "byparasp", "by10flp", "bystexp") #list of all the continuous variables

build_X <- function(df) {
  X_cat_df <- as.data.frame(df[, dummy_cols])
  X_cont   <- df |> select(all_of(vl_cont))
  cbind(X_cat_df, X_cont)
}
#^define a function to select teh continuous and dummy columns and bind them 

X_by <- build_X(df_by) #build X covariate matrix 
X_f1 <- build_X(df_f1)

cat("Feature matrix X (BY):", nrow(X_by), "rows x", ncol(X_by), "columns\n") #display observations and number of variables
cat("Feature matrix X (F1):", nrow(X_f1), "rows x", ncol(X_f1), "columns\n\n")

options(survey.lonely.psu = "adjust") #set the survey option for when ther is a single observation in a psu

dark2 <- RColorBrewer::brewer.pal(8, "Set2")
options(ggplot2.discrete.fill   = dark2,
        ggplot2.discrete.colour = dark2) #set color scheme

# Shared y range for every income CATE point-range figure, so panels are comparable
inc_ylim    <- c(-4000, 10000)
inc_ybreaks <- seq(-4000, 10000, by = 2000)

# ---- Helper: covariate balance row ----
# Computes survey-weighted means, SEs, and SMD for one variable,
# under both baseline (survey weights) and IPW-adjusted designs.
make_balance_row <- function(v, is_binary, svy_base, svy_ipw) {
  fml    <- as.formula(paste0("~`", v, "`"))
  se_col <- "se"

  base_by <- svyby(fml, ~W, svy_base, svymean, na.rm = TRUE) #survey meams by treeatment with survey weight 
  ipw_by  <- svyby(fml, ~W, svy_ipw,  svymean, na.rm = TRUE) #survey means by treatment with the combined survey*ipw weight

  w_base <- as.numeric(as.character(base_by$W))
  w_ipw  <- as.numeric(as.character(ipw_by$W))

  m0_base  <- base_by[w_base == 0, v];  se0_base <- base_by[w_base == 0, se_col]  #control
  m1_base  <- base_by[w_base == 1, v];  se1_base <- base_by[w_base == 1, se_col] #treat
  m0_ipw   <- ipw_by[w_ipw  == 0, v];  se0_ipw  <- ipw_by[w_ipw  == 0, se_col]
  m1_ipw   <- ipw_by[w_ipw  == 1, v];  se1_ipw  <- ipw_by[w_ipw  == 1, se_col]

  overall_svy <- svymean(fml, svy_base, na.rm = TRUE)
  m_overall   <- as.numeric(overall_svy)
  se_overall  <- as.numeric(SE(overall_svy)) #overall means and standard errors

  #if it is a binary variable the denominator for the smd is one formula vs the other 
  if (is_binary) {
    denom <- sqrt(max(m_overall * (1 - m_overall), 1e-10))
  } else {
    var_by <- svyby(fml, ~W, svy_base, svyvar, na.rm = TRUE)
    w_var  <- as.numeric(as.character(var_by$W))
    denom  <- sqrt(max((var_by[w_var == 0, v] + var_by[w_var == 1, v]) / 2, 1e-10))
  }

  data.frame(
    variable        = v,
    label           = get_label(v),
    mean_overall    = round(m_overall,  4),
    se_overall      = round(se_overall, 4),
    mean_treat_base = round(m1_base,    4),
    se_treat_base   = round(se1_base,   4),
    mean_ctrl_base  = round(m0_base,    4),
    se_ctrl_base    = round(se0_base,   4),
    smd_base        = round((m1_base - m0_base) / denom, 4),
    mean_treat_ipw  = round(m1_ipw,     4),
    se_treat_ipw    = round(se1_ipw,    4),
    mean_ctrl_ipw   = round(m0_ipw,     4),
    se_ctrl_ipw     = round(se0_ipw,    4),
    smd_ipw         = round((m1_ipw - m0_ipw) / denom, 4),
    row.names       = NULL
  )
}

# ---- Helper: fit one MACF + save CATE plots, return ate and forest ----
fit_macf_and_save <- function(Y_, W_, wt_, W.hat_, Y.hat_, X_macf_,
                               outdir, out_var, out_label,
                               trt_label, file_suffix, plot_suffix) {
  cat("  Fitting multi-arm causal forest (3,000 trees) [", plot_suffix, "]...\n")
  macf_ <- multi_arm_causal_forest(as.matrix(X_macf_), Y_, W_,
                                    sample.weights = wt_,
                                    W.hat = W.hat_, Y.hat = Y.hat_,
                                    num.trees = 3000)
  ate_arms_ <- average_treatment_effect(macf_, target.sample = "all")
  cat("  ATE (low vs none):",  round(ate_arms_[1, "estimate"], 5),
      "(SE:", round(ate_arms_[1, "std.err"], 5), ")\n")
  cat("  ATE (high vs none):", round(ate_arms_[2, "estimate"], 5),
      "(SE:", round(ate_arms_[2, "std.err"], 5), ")\n\n")

  tau_macf_        <- drop(predict(macf_)$predictions) #predictions return originally in the wrong number of dimensions
  contrast_labels_ <- c("Low vs None", "High vs None")
  arm_suffixes_    <- c("arm1", "arm2")
  list(ate = ate_arms_, macf = macf_, tau = tau_macf_)
}

# ---- Helper: CATE mean and 95%CI by subgroup, one panel per top-5 moderator (MACF) ----
# DR_scores_: n x 2 AIPW score matrix from get_scores(); wt_: sample weights
plot_macf_het <- function(macf_, DR_scores_, wt_, X_macf_, outdir,
                           out_var, out_label_str, trt_label, file_suffix) {
  varimp_f     <- variable_importance(macf_)
  n_top        <- min(5L, ncol(X_macf_))
  top_vars     <- colnames(X_macf_)[order(varimp_f, decreasing = TRUE)[seq_len(n_top)]]
  arm_levels   <- c("Low vs None", "High vs None")
  cat("  Het plots - top vars:", paste(top_vars, collapse = ", "), "\n")

  # Compute weighted mean and 95% CI for DR_scores[,k] within each subgroup
  group_ci <- function(dr_col, grp_vec, grp_levels) {
    do.call(rbind, lapply(grp_levels, function(g) {
      idx <- grp_vec == g & !is.na(grp_vec)
      fit <- lm(dr_col[idx] ~ 1, weights = wt_[idx])
      est <- coef(fit)[1L]
      se  <- sqrt(vcov(fit)[1L, 1L])
      data.frame(group = g, estimate = est, lo = est - 1.96 * se, hi = est + 1.96 * se)
    }))
  }

  for (v in top_vars) {
    x_v    <- X_macf_[[v]]
    is_bin <- v %in% dummy_cols || all(x_v %in% c(0, 1, NA))

    if (is_bin) {
      grp_char <- ifelse(is.na(x_v), NA_character_, ifelse(x_v == 1, "1", "0"))
      sub_labs <- c("0", "1")
    } else {
      qs       <- quantile(x_v, probs = c(0.25, 0.75), na.rm = TRUE)
      grp_char <- ifelse(x_v <= qs[1], "Quartile 1",
                  ifelse(x_v >= qs[2], "Quartile 4", NA_character_))
      sub_labs <- c("Quartile 1", "Quartile 4")
    }

    plot_df <- do.call(rbind, lapply(seq_along(arm_levels), function(k) {
      ci_df      <- group_ci(DR_scores_[, k], grp_char, sub_labs)
      ci_df$arm  <- arm_levels[k]
      ci_df
    }))
    plot_df$group <- factor(plot_df$group, levels = sub_labs)
    plot_df$arm   <- factor(plot_df$arm,   levels = arm_levels)

    grp_colors <- setNames(dark2[2:(1 + length(sub_labs))], sub_labs)
    grp_shapes <- setNames(c(16L, 17L)[seq_along(sub_labs)], sub_labs)

    p <- ggplot(plot_df, aes(x = arm, y = estimate, colour = group, shape = group)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
      geom_pointrange(aes(ymin = lo, ymax = hi),
                      position = position_dodge(width = 0.5),
                      linewidth = 0.7, size = 0.6) +
      scale_colour_manual(values = grp_colors, name = NULL) +
      scale_shape_manual(values  = grp_shapes, name = NULL) +
      labs(title    = paste("CATE by", get_label(v)),
           subtitle = paste0(trt_label, " - ", out_label_str),
           x        = "Treatment Arm",
           y        = "Average Treatment Effect") +
      theme_minimal(base_size = 11) +
      theme(legend.position    = "bottom",
            strip.text         = element_text(face = "bold"),
            panel.grid.major.x = element_blank())
    # Set consistent scale
    if (out_var == "Y_bachelors_plus") {
      p <- p + expand_limits(y = c(-0.25, 0.15)) +
        scale_y_continuous(breaks = seq(-0.2, 0.2, by = 0.1))
    } else {
      p <- p + scale_y_continuous(limits = inc_ylim, breaks = inc_ybreaks,
                                  labels = scales::label_comma())
      if (min(plot_df$lo) < inc_ylim[1] || max(plot_df$hi) > inc_ylim[2]) {
        cat("    WARNING: CI outside inc_ylim for", v,
            "- observed [", round(min(plot_df$lo)), ",",
            round(max(plot_df$hi)), "]\n")
      }
    }
    
   
    v_safe <- gsub("[^A-Za-z0-9_]", "_", v)
    fname  <- paste0("cate_het_", out_var, "_", v_safe, file_suffix, ".tiff")
    ggsave(file.path(outdir, fname), p, width = 6, height = 4, dpi = 600)
    cat("    Saved:", fname, "\n")
  }
}

# ==============================================================================
# 1. BY treatment — any work in 10th grade (D_by)
# ==============================================================================

cat("=====================================\n")
cat("Treatment: Work (10th grade, D_by)\n")
cat("=====================================\n\n")

#Follows a general structure: fit propensity score for both outcomes; plot overlap of pscore
#create balance tables; fit outcome forest for BA; fit initial CF for variable selection
# fit final CF on top 1/3 of variables; save CATE; 50/50 split for RATE estimation
# run Best linear prediction IF RATE is significant. Then repeat for other outcome. 

outdir_by <- "output/by_cf" #set output directory

W_by  <- as.numeric(df_by$D_by) #Treatment variable
wt_by <- df_by$f3bypnlwt        #Set the weight
cat("N:", nrow(df_by), "| Treated:", sum(W_by), "| Control:", sum(W_by == 0), "\n\n") #display

#  Propensity forest (shared across outcomes) -----------------------------
cat("Fitting propensity forest (BY)...\n")
forest.W_by <- probability_forest(X_by, factor(W_by), sample.weights = wt_by,
                                   num.trees = 2000, honesty = TRUE) #only accepts factor variables for treatment
W.hat_by <- predict(forest.W_by)$predictions[, "1"] #store the predicted treatment probability (propensity score)


cat("W.hat_by range: [", round(min(W.hat_by), 3), ",", round(max(W.hat_by), 3), "]\n\n") #display 10th grade propensity range

# Overlap plot
p_overlap_by <- ggplot(data.frame(W.hat = W.hat_by,
                                   D = factor(W_by, levels = c(0, 1),
                                              labels = c("No work", "Work"))),
                        aes(x = W.hat, fill = D, color = D)) +
  geom_density(alpha = 0.6) +
  labs(fill = NULL, color = NULL, x = "Estimated Propensity Score",
       title = "Overlap: Work (10th grade)") +
    expand_limits(x = c(0, 1)) + #set axis so it has to span 0 to 1
  theme_minimal()
ggsave(file.path(outdir_by, "overlap.tiff"), p_overlap_by, width = 6, height = 4, dpi = 600)

# IPW and combined weight
IPW_by    <- ifelse(W_by == 1, 1 / W.hat_by, 1 / (1 - W.hat_by)) #IPW weight is 1/propensity if treated, 1-phat if untreated
w_comb_by <- wt_by * IPW_by

# Survey designs for balance table
X_cat_by <- X_by[, dummy_cols] #all the cateforical variable 

svy_base_by <- svydesign(
  id = ~psu, strata = ~strat_id, nest = TRUE,
  weights = ~wt_base,
  data = cbind(df_by[, c("psu", "strat_id")],
               data.frame(W = W_by, wt_base = wt_by),
               X_cat_by, df_by[, vl_cont])
)
svy_ipw_by <- svydesign(
  id = ~psu, strata = ~strat_id, nest = TRUE,
  weights = ~w_comb,
  data = cbind(df_by[, c("psu", "strat_id")],
               data.frame(W = W_by, w_comb = w_comb_by),
               X_cat_by, df_by[, vl_cont])
)

# Balance table — lapply over covariates 
bal_bin_by  <- Filter(Negate(is.null), lapply(names(X_cat_by), function(v) {
  tryCatch(make_balance_row(v, is_binary = TRUE,  svy_base_by, svy_ipw_by), error = function(e) NULL)
})) #over the categoricals, save a list of data frames with all the stats
bal_cont_by <- Filter(Negate(is.null), lapply(vl_cont, function(v) {
  tryCatch(make_balance_row(v, is_binary = FALSE, svy_base_by, svy_ipw_by), error = function(e) NULL)
})) #Over the continuous
balance_by <- do.call(rbind, c(bal_bin_by, bal_cont_by)) #row bind aka append
write.csv(balance_by, file.path(outdir_by, "balance_table.csv"), row.names = FALSE)
cat("Balance table saved:", file.path(outdir_by, "balance_table.csv"), "\n\n")

# ---- Outcome: BA attainment (BY)-----------------------
cat("--- Outcome: BA or Higher (BY) ---\n")
Y_by_bach <- as.numeric(df_by$Y_bachelors_plus) #bachelors outcome

cat("  Fitting outcome forest (BY, bachelors)...\n")
forest.Y_by_bach <- regression_forest(X_by, Y_by_bach, sample.weights = wt_by,
                                       num.trees = 1500, tune.parameters = "all")
#predict the outcome 
Y.hat_by_bach <- predict(forest.Y_by_bach)$predictions #store the predicted values

cat("  Fitting initial causal forest for variable selection (BY, bachelors)...\n")
forest.init_by_bach <- causal_forest(as.matrix(X_by), Y_by_bach, W_by,
                                      sample.weights = wt_by,
                                      Y.hat = Y.hat_by_bach, W.hat = W.hat_by,
                                      num.trees = 2000, tune.parameters = "none")
varimp_by_bach <- variable_importance(forest.init_by_bach)

n_keep_by_bach <- max(1L, round(ncol(X_by) / 5)) #keep the top fifth of most important variable 
keep_by_bach   <- colnames(X_by)[order(varimp_by_bach, decreasing = TRUE)[seq_len(n_keep_by_bach)]]
cat("  Top fifth vars (", n_keep_by_bach, "):", paste(keep_by_bach, collapse = ", "), "\n") #print them in the log

X_cate_by_bach <- X_by[, keep_by_bach]
cat("  Fitting causal forest (6,000 trees) (BY, bachelors)...\n")
cf_by_bach <- causal_forest(as.matrix(X_cate_by_bach), Y_by_bach, W_by,
                             sample.weights = wt_by,
                             Y.hat = Y.hat_by_bach, W.hat = W.hat_by,
                             num.trees = 6000, tune.parameters = "all")

varimp_final_by_bach <- variable_importance(cf_by_bach) #store final variable import
cf_pred_by_bach <- predict(cf_by_bach) #get the predictions

ate_by_bach <- average_treatment_effect(cf_by_bach, target.sample = "all")
att_by_bach <- average_treatment_effect(cf_by_bach, target.sample = "treated")
cat("  ATE:", round(ate_by_bach["estimate"], 5), "(SE:", round(ate_by_bach["std.err"], 5), ")\n")
cat("  ATT:", round(att_by_bach["estimate"], 5), "(SE:", round(att_by_bach["std.err"], 5), ")\n")

cat("  Calibration test:\n")
cal_by_bach <- test_calibration(cf_by_bach)
print(cal_by_bach)

# Rank average treatment effect -  50/50 sample split
cat("  Computing RATE via 50/50 sample splitting (BY, bachelors)...\n")
set.seed(1)
n_by_bach      <- nrow(X_cate_by_bach) #the the number of observations

train_idx_by_bach <- sample(seq_len(n_by_bach), floor(n_by_bach / 2)) #Sample the data frame
eval_idx_by_bach  <- setdiff(seq_len(n_by_bach), train_idx_by_bach) #keep the other half
cat("  Split: train =", length(train_idx_by_bach), "| eval =", length(eval_idx_by_bach), "\n")

p_by_bach <- cf_by_bach$tunable.params #pass the tuning parameters to save time
cf_train_by_bach <- causal_forest(
  as.matrix(X_cate_by_bach[train_idx_by_bach, ]), Y_by_bach[train_idx_by_bach], W_by[train_idx_by_bach],
  sample.weights = wt_by[train_idx_by_bach], num.trees = 3000, tune.parameters = "none",
  min.node.size = p_by_bach$min.node.size, sample.fraction = p_by_bach$sample.fraction,
  mtry = p_by_bach$mtry, alpha = p_by_bach$alpha, imbalance.penalty = p_by_bach$imbalance.penalty)
cf_eval_by_bach  <- causal_forest(
  as.matrix(X_cate_by_bach[eval_idx_by_bach, ]),  Y_by_bach[eval_idx_by_bach],  W_by[eval_idx_by_bach],
  sample.weights = wt_by[eval_idx_by_bach],  num.trees = 3000, tune.parameters = "none",
  min.node.size = p_by_bach$min.node.size, sample.fraction = p_by_bach$sample.fraction,
  mtry = p_by_bach$mtry, alpha = p_by_bach$alpha, imbalance.penalty = p_by_bach$imbalance.penalty)
priorities_by_bach <- predict(cf_train_by_bach,
                               as.matrix(X_cate_by_bach[eval_idx_by_bach, ]))$predictions

rate_autoc_by_bach <- rank_average_treatment_effect(cf_eval_by_bach, priorities_by_bach, target = "AUTOC")
rate_autoc_p_by_bach <- 2 * pnorm(-abs(rate_autoc_by_bach$estimate / rate_autoc_by_bach$std.err))

cat("  RATE AUTOC:", round(rate_autoc_by_bach$estimate, 5),
    "(SE:", round(rate_autoc_by_bach$std.err, 5), "| p:", round(rate_autoc_p_by_bach, 4), ")\n")

toc_path_by_bach <- file.path(outdir_by, "rate_toc_Y_bachelors_plus.tiff")
tiff(toc_path_by_bach, width = 6, height = 4, units = "in", res = 600)
plot(rate_autoc_by_bach,
     main = "TOC Curve (AUTOC, held-out half): Work (10th grade)\nBA or Higher")
dev.off()

rate_sig_by_bach <- (rate_autoc_p_by_bach < 0.05) #This is a gate, doing BLP only if significant
blp_df_by_bach <- NULL
if (rate_sig_by_bach) {
  cat("  RATE significant: fitting BLP on top 5 variables (BY, bachelors)...\n")

  n_blp_by_bach    <- min(5L, ncol(X_cate_by_bach))
  blp_vars_by_bach <- colnames(X_cate_by_bach)[
    order(varimp_final_by_bach, decreasing = TRUE)[seq_len(n_blp_by_bach)]
  ]
  cat("  BLP vars:", paste(blp_vars_by_bach, collapse = ", "), "\n")

  blp_df_by_bach <- do.call(rbind, lapply(blp_vars_by_bach, function(var) {
    res <- best_linear_projection(cf_by_bach,
                                  A = as.matrix(X_cate_by_bach[, var, drop = FALSE]))
    print(res)
    m <- as.data.frame(unclass(res))
    data.frame(source_variable  = var,
               term             = rownames(m),
               Estimate         = round(m$Estimate,         6),
               Std_Error        = round(m[["Std. Error"]],  6),
               t_value          = round(m[["t value"]],     4),
               p_value          = m[["Pr(>|t|)"]],
               stringsAsFactors = FALSE, row.names = NULL)
  }))
  slope_rows <- blp_df_by_bach$term != "(Intercept)"
  blp_df_by_bach$p_adj_BH <- NA_real_
  blp_df_by_bach$p_adj_BH[slope_rows] <- p.adjust(blp_df_by_bach$p_value[slope_rows], method = "BH")
  write.csv(blp_df_by_bach, file = file.path(outdir_by, "blp_BA.csv"), row.names = FALSE)
  cat("  BLP saved:", file.path(outdir_by, "blp_BA.csv"), "\n")
}

#  Outcome: Annual earnings (BY)
cat("--- Outcome: Annual Earnings (BY) ---\n")
Y_by_inc <- as.numeric(df_by$Y_income)

cat("  Fitting outcome forest (BY, income)...\n")
forest.Y_by_inc <- regression_forest(X_by, Y_by_inc, sample.weights = wt_by,
                                      num.trees = 1500, tune.parameters = "all") #train outcome forest
Y.hat_by_inc <- predict(forest.Y_by_inc)$predictions

cat("  Fitting initial causal forest for variable selection (BY, income)...\n")
forest.init_by_inc <- causal_forest(as.matrix(X_by), Y_by_inc, W_by,
                                     sample.weights = wt_by,
                                     Y.hat = Y.hat_by_inc, W.hat = W.hat_by,
                                     num.trees = 2000, tune.parameters = "none")
varimp_by_inc <- variable_importance(forest.init_by_inc)
n_keep_by_inc <- max(1L, round(ncol(X_by) / 5))
keep_by_inc   <- colnames(X_by)[order(varimp_by_inc, decreasing = TRUE)[seq_len(n_keep_by_inc)]]
cat("  Top fifth vars (", n_keep_by_inc, "):", paste(keep_by_inc, collapse = ", "), "\n")

X_cate_by_inc <- X_by[, keep_by_inc]
cat("  Fitting causal forest (6,000 trees) (BY, income)...\n")
cf_by_inc <- causal_forest(as.matrix(X_cate_by_inc), Y_by_inc, W_by,
                            sample.weights = wt_by,
                            Y.hat = Y.hat_by_inc, W.hat = W.hat_by,
                            num.trees = 6000, tune.parameters = "all")
varimp_final_by_inc <- variable_importance(cf_by_inc)
cf_pred_by_inc <- predict(cf_by_inc)

ate_by_inc <- average_treatment_effect(cf_by_inc, target.sample = "all")
att_by_inc <- average_treatment_effect(cf_by_inc, target.sample = "treated")
cat("  ATE:", round(ate_by_inc["estimate"], 5), "(SE:", round(ate_by_inc["std.err"], 5), ")\n")
cat("  ATT:", round(att_by_inc["estimate"], 5), "(SE:", round(att_by_inc["std.err"], 5), ")\n")

cat("  Computing RATE via 50/50 sample splitting (BY, income)...\n")
set.seed(1)
n_by_inc         <- nrow(X_cate_by_inc)
train_idx_by_inc <- sample(seq_len(n_by_inc), floor(n_by_inc / 2))
eval_idx_by_inc  <- setdiff(seq_len(n_by_inc), train_idx_by_inc)
cat("  Split: train =", length(train_idx_by_inc), "| eval =", length(eval_idx_by_inc), "\n")

p_by_inc <- cf_by_inc$tunable.params
cf_train_by_inc <- causal_forest(
  as.matrix(X_cate_by_inc[train_idx_by_inc, ]), Y_by_inc[train_idx_by_inc], W_by[train_idx_by_inc],
  sample.weights = wt_by[train_idx_by_inc], num.trees = 3000, tune.parameters = "none",
  min.node.size = p_by_inc$min.node.size, sample.fraction = p_by_inc$sample.fraction,
  mtry = p_by_inc$mtry, alpha = p_by_inc$alpha, imbalance.penalty = p_by_inc$imbalance.penalty)
cf_eval_by_inc  <- causal_forest(
  as.matrix(X_cate_by_inc[eval_idx_by_inc, ]),  Y_by_inc[eval_idx_by_inc],  W_by[eval_idx_by_inc],
  sample.weights = wt_by[eval_idx_by_inc],  num.trees = 3000, tune.parameters = "none",
  min.node.size = p_by_inc$min.node.size, sample.fraction = p_by_inc$sample.fraction,
  mtry = p_by_inc$mtry, alpha = p_by_inc$alpha, imbalance.penalty = p_by_inc$imbalance.penalty)
priorities_by_inc <- predict(cf_train_by_inc,
                              as.matrix(X_cate_by_inc[eval_idx_by_inc, ]))$predictions

rate_autoc_by_inc <- rank_average_treatment_effect(cf_eval_by_inc, priorities_by_inc, target = "AUTOC")

rate_autoc_p_by_inc <- 2 * pnorm(-abs(rate_autoc_by_inc$estimate / rate_autoc_by_inc$std.err))

cat("  RATE AUTOC:", round(rate_autoc_by_inc$estimate, 5),
    "(SE:", round(rate_autoc_by_inc$std.err, 5), "| p:", round(rate_autoc_p_by_inc, 4), ")\n")

toc_path_by_inc <- file.path(outdir_by, "rate_toc_Y_income.tiff")
tiff(toc_path_by_inc, width = 6, height = 4, units = "in", res = 600)
plot(rate_autoc_by_inc,
     main = "TOC Curve (AUTOC, held-out half): Work (10th grade)\nAnnual Earnings ($)")
dev.off()

rate_sig_by_inc <- (rate_autoc_p_by_inc < 0.05)
blp_df_by_inc <- NULL
if (rate_sig_by_inc) {
  cat("  RATE significant: fitting BLP on top 5 variables (BY, income)...\n")

  n_blp_by_inc    <- min(5L, ncol(X_cate_by_inc))
  blp_vars_by_inc <- colnames(X_cate_by_inc)[
    order(varimp_final_by_inc, decreasing = TRUE)[seq_len(n_blp_by_inc)]
  ]
  cat("  BLP vars:", paste(blp_vars_by_inc, collapse = ", "), "\n")

  blp_df_by_inc <- do.call(rbind, lapply(blp_vars_by_inc, function(var) {
    res <- best_linear_projection(cf_by_inc,
                                  A = as.matrix(X_cate_by_inc[, var, drop = FALSE]))
    print(res)
    m <- as.data.frame(unclass(res))
    data.frame(source_variable  = var,
               term             = rownames(m),
               Estimate         = round(m$Estimate,         6),
               Std_Error        = round(m[["Std. Error"]],  6),
               t_value          = round(m[["t value"]],     4),
               p_value          = m[["Pr(>|t|)"]],
               stringsAsFactors = FALSE, row.names = NULL)
  }))
  slope_rows <- blp_df_by_inc$term != "(Intercept)"
  blp_df_by_inc$p_adj_BH <- NA_real_
  blp_df_by_inc$p_adj_BH[slope_rows] <- p.adjust(blp_df_by_inc$p_value[slope_rows], method = "BH")
  write.csv(blp_df_by_inc, file = file.path(outdir_by, "blp_INC.csv"), row.names = FALSE)
  cat("  BLP saved:", file.path(outdir_by, "blp_INC.csv"), "\n")
}


# -- Results table: BY treatment -------------------------
# ATE / ATT
ate_rows_by <- data.frame(
  Outcome = c("BA or Higher", "Annual Earnings ($)"),
  ATE     = round(c(ate_by_bach["estimate"], ate_by_inc["estimate"]), 4),
  ATE_SE  = round(c(ate_by_bach["std.err"],  ate_by_inc["std.err"]),  4),
  ATT     = round(c(att_by_bach["estimate"], att_by_inc["estimate"]), 4),
  ATT_SE  = round(c(att_by_bach["std.err"],  att_by_inc["std.err"]),  4),
  row.names = NULL
)

#put results in csv 
write.csv(ate_rows_by, file = file.path(outdir_by,"results.csv"), row.names = TRUE)

# RATE
rate_rows_by <- data.frame(
  Outcome  = c("BA or Higher", "Annual Earnings ($)"),
  Metric   = c("AUTOC",  "AUTOC"),
  Estimate = round(c(rate_autoc_by_bach$estimate, 
                     rate_autoc_by_inc$estimate), 5),
  SE       = round(c(rate_autoc_by_bach$std.err,  
                     rate_autoc_by_inc$std.err),  5),
  p_value  = round(c(rate_autoc_p_by_bach,
                     rate_autoc_p_by_inc), 4),
  row.names = NULL
)
#put results in csv 
write.csv(rate_rows_by, file = file.path(outdir_by,"rate.csv"), row.names = TRUE)


cat("Results table saved, rate table, blp if significant", "\n\n")

# ==============================================================================
# 2. F1 treatment — any work in 12th grade (D_f1)
# ==============================================================================

cat("=====================================\n")
cat("Treatment: Work (12th grade, D_f1)\n")
cat("=====================================\n\n")

outdir_f1 <- "output/f1_cf"

W_f1  <- as.numeric(df_f1$D_f1)
wt_f1 <- df_f1$f3bypnlwt
cat("N:", nrow(df_f1), "| Treated:", sum(W_f1), "| Control:", sum(W_f1 == 0), "\n\n")

#  Propensity forest (shared across outcomes) 
cat("Fitting propensity forest (F1)...\n")
forest.W_f1 <- probability_forest(X_f1, factor(W_f1), sample.weights = wt_f1,
                                   num.trees = 2000, honesty = TRUE)
W.hat_f1 <- predict(forest.W_f1)$predictions[, "1"]
cat("W.hat_f1 range: [", round(min(W.hat_f1), 3), ",", round(max(W.hat_f1), 3), "]\n\n")

p_overlap_f1 <- ggplot(data.frame(W.hat = W.hat_f1,
                                   D = factor(W_f1, levels = c(0, 1),
                                              labels = c("No work", "Work"))),
                        aes(x = W.hat, fill = D, color = D)) +
  geom_density(alpha = 0.6) +
  labs(fill = NULL, color = NULL, x = "Estimated Propensity Score",
       title = "Overlap: Work (12th grade)") +
  expand_limits(x = c(0, 1)) + 
  theme_minimal()
ggsave(file.path(outdir_f1, "overlap.tiff"), p_overlap_f1, width = 6, height = 4, dpi = 600)

IPW_f1    <- ifelse(W_f1 == 1, 1 / W.hat_f1, 1 / (1 - W.hat_f1))
w_comb_f1 <- wt_f1 * IPW_f1

X_cat_f1 <- X_f1[, dummy_cols]

svy_base_f1 <- svydesign(
  id = ~psu, strata = ~strat_id, nest = TRUE,
  weights = ~wt_base,
  data = cbind(df_f1[, c("psu", "strat_id")],
               data.frame(W = W_f1, wt_base = wt_f1),
               X_cat_f1, df_f1[, vl_cont])
)
svy_ipw_f1 <- svydesign(
  id = ~psu, strata = ~strat_id, nest = TRUE,
  weights = ~w_comb,
  data = cbind(df_f1[, c("psu", "strat_id")],
               data.frame(W = W_f1, w_comb = w_comb_f1),
               X_cat_f1, df_f1[, vl_cont])
)

bal_bin_f1  <- Filter(Negate(is.null), lapply(names(X_cat_f1), function(v) {
  tryCatch(make_balance_row(v, is_binary = TRUE,  svy_base_f1, svy_ipw_f1), error = function(e) NULL)
}))
bal_cont_f1 <- Filter(Negate(is.null), lapply(vl_cont, function(v) {
  tryCatch(make_balance_row(v, is_binary = FALSE, svy_base_f1, svy_ipw_f1), error = function(e) NULL)
}))
balance_f1 <- do.call(rbind, c(bal_bin_f1, bal_cont_f1))
write.csv(balance_f1, file.path(outdir_f1, "balance_table.csv"), row.names = FALSE)
cat("Balance table saved:", file.path(outdir_f1, "balance_table.csv"), "\n\n")

# ── Outcome: BA attainment (F1) ──────────────────────────────────────────────
cat("--- Outcome: BA or Higher (F1) ---\n")
Y_f1_bach <- as.numeric(df_f1$Y_bachelors_plus)

cat("  Fitting outcome forest (F1, bachelors)...\n")
forest.Y_f1_bach <- regression_forest(X_f1, Y_f1_bach, sample.weights = wt_f1,
                                       num.trees = 1500, tune.parameters = "all")
Y.hat_f1_bach <- predict(forest.Y_f1_bach)$predictions

cat("  Fitting initial causal forest for variable selection (F1, bachelors)...\n")
forest.init_f1_bach <- causal_forest(as.matrix(X_f1), Y_f1_bach, W_f1,
                                      sample.weights = wt_f1,
                                      Y.hat = Y.hat_f1_bach, W.hat = W.hat_f1,
                                      num.trees = 2000, tune.parameters = "none")
varimp_f1_bach <- variable_importance(forest.init_f1_bach)
n_keep_f1_bach <- max(1L, round(ncol(X_f1) / 5))
keep_f1_bach   <- colnames(X_f1)[order(varimp_f1_bach, decreasing = TRUE)[seq_len(n_keep_f1_bach)]]
cat("  Top fifth vars (", n_keep_f1_bach, "):", paste(keep_f1_bach, collapse = ", "), "\n")

X_cate_f1_bach <- X_f1[, keep_f1_bach]
cat("  Fitting causal forest (6,000 trees) (F1, bachelors)...\n")
cf_f1_bach <- causal_forest(as.matrix(X_cate_f1_bach), Y_f1_bach, W_f1,
                             sample.weights = wt_f1,
                             Y.hat = Y.hat_f1_bach, W.hat = W.hat_f1,
                             num.trees = 6000, tune.parameters = "all")
varimp_final_f1_bach <- variable_importance(cf_f1_bach)
cf_pred_f1_bach <- predict(cf_f1_bach)

ate_f1_bach <- average_treatment_effect(cf_f1_bach, target.sample = "all")
att_f1_bach <- average_treatment_effect(cf_f1_bach, target.sample = "treated")
cat("  ATE:", round(ate_f1_bach["estimate"], 5), "(SE:", round(ate_f1_bach["std.err"], 5), ")\n")
cat("  ATT:", round(att_f1_bach["estimate"], 5), "(SE:", round(att_f1_bach["std.err"], 5), ")\n")

cat("  Computing RATE via 50/50 sample splitting (F1, bachelors)...\n")
set.seed(1)
n_f1_bach         <- nrow(X_cate_f1_bach)
train_idx_f1_bach <- sample(seq_len(n_f1_bach), floor(n_f1_bach / 2))
eval_idx_f1_bach  <- setdiff(seq_len(n_f1_bach), train_idx_f1_bach)
cat("  Split: train =", length(train_idx_f1_bach), "| eval =", length(eval_idx_f1_bach), "\n")

p_f1_bach <- cf_f1_bach$tunable.params
cf_train_f1_bach <- causal_forest(
  as.matrix(X_cate_f1_bach[train_idx_f1_bach, ]), Y_f1_bach[train_idx_f1_bach], W_f1[train_idx_f1_bach],
  sample.weights = wt_f1[train_idx_f1_bach], num.trees = 3000, tune.parameters = "none",
  min.node.size = p_f1_bach$min.node.size, sample.fraction = p_f1_bach$sample.fraction,
  mtry = p_f1_bach$mtry, alpha = p_f1_bach$alpha, imbalance.penalty = p_f1_bach$imbalance.penalty)
cf_eval_f1_bach  <- causal_forest(
  as.matrix(X_cate_f1_bach[eval_idx_f1_bach, ]),  Y_f1_bach[eval_idx_f1_bach],  W_f1[eval_idx_f1_bach],
  sample.weights = wt_f1[eval_idx_f1_bach],  num.trees = 3000, tune.parameters = "none",
  min.node.size = p_f1_bach$min.node.size, sample.fraction = p_f1_bach$sample.fraction,
  mtry = p_f1_bach$mtry, alpha = p_f1_bach$alpha, imbalance.penalty = p_f1_bach$imbalance.penalty)
priorities_f1_bach <- predict(cf_train_f1_bach,
                               as.matrix(X_cate_f1_bach[eval_idx_f1_bach, ]))$predictions

rate_autoc_f1_bach <- rank_average_treatment_effect(cf_eval_f1_bach, priorities_f1_bach, target = "AUTOC")
rate_autoc_p_f1_bach <- 2 * pnorm(-abs(rate_autoc_f1_bach$estimate / rate_autoc_f1_bach$std.err))

cat("  RATE AUTOC:", round(rate_autoc_f1_bach$estimate, 5),
    "(SE:", round(rate_autoc_f1_bach$std.err, 5), "| p:", round(rate_autoc_p_f1_bach, 4), ")\n")

toc_path_f1_bach <- file.path(outdir_f1, "rate_toc_Y_bachelors_plus.tiff")
tiff(toc_path_f1_bach, width = 6, height = 4, units = "in", res = 600)
plot(rate_autoc_f1_bach,
     main = "TOC Curve (AUTOC, held-out half): Work (12th grade)\nBA or Higher")
dev.off()

rate_sig_f1_bach <- (rate_autoc_p_f1_bach < 0.05)
blp_df_f1_bach <- NULL
if (rate_sig_f1_bach) {
  cat("  RATE significant: fitting BLP on top 5 variables (F1, bachelors)...\n")

  n_blp_f1_bach    <- min(5L, ncol(X_cate_f1_bach))
  blp_vars_f1_bach <- colnames(X_cate_f1_bach)[
    order(varimp_final_f1_bach, decreasing = TRUE)[seq_len(n_blp_f1_bach)]
  ]
  cat("  BLP vars:", paste(blp_vars_f1_bach, collapse = ", "), "\n")

  blp_df_f1_bach <- do.call(rbind, lapply(blp_vars_f1_bach, function(var) {
    res <- best_linear_projection(cf_f1_bach,
                                  A = as.matrix(X_cate_f1_bach[, var, drop = FALSE]))
    print(res)
    m <- as.data.frame(unclass(res))
    data.frame(source_variable  = var,
               term             = rownames(m),
               Estimate         = round(m$Estimate,         6),
               Std_Error        = round(m[["Std. Error"]],  6),
               t_value          = round(m[["t value"]],     4),
               p_value          = m[["Pr(>|t|)"]],
               stringsAsFactors = FALSE, row.names = NULL)
  }))
  slope_rows <- blp_df_f1_bach$term != "(Intercept)"
  blp_df_f1_bach$p_adj_BH <- NA_real_
  blp_df_f1_bach$p_adj_BH[slope_rows] <- p.adjust(blp_df_f1_bach$p_value[slope_rows], method = "BH")
  write.csv(blp_df_f1_bach, file = file.path(outdir_f1, "blp_BA.csv"), row.names = FALSE)
  cat("  BLP saved:", file.path(outdir_f1, "blp_BA.csv"), "\n")
}

# ── Outcome: Annual earnings (F1) ────────────────────────────────────────────
cat("--- Outcome: Annual Earnings (F1) ---\n")
Y_f1_inc <- as.numeric(df_f1$Y_income)

cat("  Fitting outcome forest (F1, income)...\n")
forest.Y_f1_inc <- regression_forest(X_f1, Y_f1_inc, sample.weights = wt_f1,
                                      num.trees = 1500, tune.parameters = "all")
Y.hat_f1_inc <- predict(forest.Y_f1_inc)$predictions

cat("  Fitting initial causal forest for variable selection (F1, income)...\n")
forest.init_f1_inc <- causal_forest(as.matrix(X_f1), Y_f1_inc, W_f1,
                                     sample.weights = wt_f1,
                                     Y.hat = Y.hat_f1_inc, W.hat = W.hat_f1,
                                     num.trees = 2000, tune.parameters = "none")
varimp_f1_inc <- variable_importance(forest.init_f1_inc)
n_keep_f1_inc <- max(1L, round(ncol(X_f1) / 5))
keep_f1_inc   <- colnames(X_f1)[order(varimp_f1_inc, decreasing = TRUE)[seq_len(n_keep_f1_inc)]]
cat("  Top fifth vars (", n_keep_f1_inc, "):", paste(keep_f1_inc, collapse = ", "), "\n")

X_cate_f1_inc <- X_f1[, keep_f1_inc]
cat("  Fitting causal forest (6,000 trees) (F1, income)...\n")
cf_f1_inc <- causal_forest(as.matrix(X_cate_f1_inc), Y_f1_inc, W_f1,
                            sample.weights = wt_f1,
                            Y.hat = Y.hat_f1_inc, W.hat = W.hat_f1,
                            num.trees = 6000, tune.parameters = "all")
varimp_final_f1_inc <- variable_importance(cf_f1_inc)
cf_pred_f1_inc <- predict(cf_f1_inc)

ate_f1_inc <- average_treatment_effect(cf_f1_inc, target.sample = "all")
att_f1_inc <- average_treatment_effect(cf_f1_inc, target.sample = "treated")
cat("  ATE:", round(ate_f1_inc["estimate"], 5), "(SE:", round(ate_f1_inc["std.err"], 5), ")\n")
cat("  ATT:", round(att_f1_inc["estimate"], 5), "(SE:", round(att_f1_inc["std.err"], 5), ")\n")

cat("  Computing RATE via 50/50 sample splitting (F1, income)...\n")
set.seed(1)
n_f1_inc         <- nrow(X_cate_f1_inc)
train_idx_f1_inc <- sample(seq_len(n_f1_inc), floor(n_f1_inc / 2))
eval_idx_f1_inc  <- setdiff(seq_len(n_f1_inc), train_idx_f1_inc)
cat("  Split: train =", length(train_idx_f1_inc), "| eval =", length(eval_idx_f1_inc), "\n")

p_f1_inc <- cf_f1_inc$tunable.params
cf_train_f1_inc <- causal_forest(
  as.matrix(X_cate_f1_inc[train_idx_f1_inc, ]), Y_f1_inc[train_idx_f1_inc], W_f1[train_idx_f1_inc],
  sample.weights = wt_f1[train_idx_f1_inc], num.trees = 3000, tune.parameters = "none",
  min.node.size = p_f1_inc$min.node.size, sample.fraction = p_f1_inc$sample.fraction,
  mtry = p_f1_inc$mtry, alpha = p_f1_inc$alpha, imbalance.penalty = p_f1_inc$imbalance.penalty)
cf_eval_f1_inc  <- causal_forest(
  as.matrix(X_cate_f1_inc[eval_idx_f1_inc, ]),  Y_f1_inc[eval_idx_f1_inc],  W_f1[eval_idx_f1_inc],
  sample.weights = wt_f1[eval_idx_f1_inc],  num.trees = 3000, tune.parameters = "none",
  min.node.size = p_f1_inc$min.node.size, sample.fraction = p_f1_inc$sample.fraction,
  mtry = p_f1_inc$mtry, alpha = p_f1_inc$alpha, imbalance.penalty = p_f1_inc$imbalance.penalty)
priorities_f1_inc <- predict(cf_train_f1_inc,
                              as.matrix(X_cate_f1_inc[eval_idx_f1_inc, ]))$predictions

rate_autoc_f1_inc <- rank_average_treatment_effect(cf_eval_f1_inc, priorities_f1_inc, target = "AUTOC")
rate_autoc_p_f1_inc <- 2 * pnorm(-abs(rate_autoc_f1_inc$estimate / rate_autoc_f1_inc$std.err))
cat("  RATE AUTOC:", round(rate_autoc_f1_inc$estimate, 5),
    "(SE:", round(rate_autoc_f1_inc$std.err, 5), "| p:", round(rate_autoc_p_f1_inc, 4), ")\n")

toc_path_f1_inc <- file.path(outdir_f1, "rate_toc_Y_income.tiff")
tiff(toc_path_f1_inc, width = 6, height = 4, units = "in", res = 600)
plot(rate_autoc_f1_inc,
     main = "TOC Curve (AUTOC, held-out half): Work (12th grade)\nAnnual Earnings ($)")
dev.off()

rate_sig_f1_inc <- (rate_autoc_p_f1_inc < 0.05)
blp_df_f1_inc <- NULL
if (rate_sig_f1_inc) {
  cat("  RATE significant: fitting BLP on top 5 variables (F1, income)...\n")

  n_blp_f1_inc    <- min(5L, ncol(X_cate_f1_inc))
  blp_vars_f1_inc <- colnames(X_cate_f1_inc)[
    order(varimp_final_f1_inc, decreasing = TRUE)[seq_len(n_blp_f1_inc)]
  ]
  cat("  BLP vars:", paste(blp_vars_f1_inc, collapse = ", "), "\n")

  blp_df_f1_inc <- do.call(rbind, lapply(blp_vars_f1_inc, function(var) {
    res <- best_linear_projection(cf_f1_inc,
                                  A = as.matrix(X_cate_f1_inc[, var, drop = FALSE]))
    print(res)
    m <- as.data.frame(unclass(res))
    data.frame(source_variable  = var,
               term             = rownames(m),
               Estimate         = round(m$Estimate,         6),
               Std_Error        = round(m[["Std. Error"]],  6),
               t_value          = round(m[["t value"]],     4),
               p_value          = m[["Pr(>|t|)"]],
               stringsAsFactors = FALSE, row.names = NULL)
  }))
  slope_rows <- blp_df_f1_inc$term != "(Intercept)"
  blp_df_f1_inc$p_adj_BH <- NA_real_
  blp_df_f1_inc$p_adj_BH[slope_rows] <- p.adjust(blp_df_f1_inc$p_value[slope_rows], method = "BH")
  write.csv(blp_df_f1_inc, file = file.path(outdir_f1, "blp_INC.csv"), row.names = FALSE)
  cat("  BLP saved:", file.path(outdir_f1, "blp_INC.csv"), "\n")
}

# ── Results table: F1 treatment ───────────────────────────────────────────────
ate_rows_f1 <- data.frame(
  Outcome = c("BA or Higher", "Annual Earnings ($)"),
  ATE     = round(c(ate_f1_bach["estimate"], ate_f1_inc["estimate"]), 4),
  ATE_SE  = round(c(ate_f1_bach["std.err"],  ate_f1_inc["std.err"]),  4),
  ATT     = round(c(att_f1_bach["estimate"], att_f1_inc["estimate"]), 4),
  ATT_SE  = round(c(att_f1_bach["std.err"],  att_f1_inc["std.err"]),  4),
  row.names = NULL
)

write.csv(ate_rows_f1, file = file.path(outdir_f1, "results.csv"), row.names = TRUE)


# RATE
rate_rows_f1 <- data.frame(
  Outcome  = c("BA or Higher", "Annual Earnings ($)"),
  Metric   = c("AUTOC",  "AUTOC"),
  Estimate = round(c(rate_autoc_f1_bach$estimate, 
                     rate_autoc_f1_inc$estimate), 5),
  SE       = round(c(rate_autoc_f1_bach$std.err,  
                     rate_autoc_f1_inc$std.err),  5),
  p_value  = round(c(rate_autoc_p_f1_bach,
                     rate_autoc_p_f1_inc), 4),
  row.names = NULL
)
#put results in csv 
write.csv(rate_rows_f1, file = file.path(outdir_f1,"rate.csv"), row.names = TRUE)


cat("Results table saved, rate table, blp if significant", "\n\n")

# ==============================================================================
# 3. MACF — work intensity in 12th grade (W_f1: 0=none, 1=1-20 hrs/wk, 2=21+)
# ==============================================================================

cat("=====================================\n")
cat("MACF: Work Intensity (12th grade, W_f1)\n")
cat("=====================================\n\n")

# Subset to non-missing intensity treatment
macf_idx  <- !is.na(df_f1$W_f1)
n_dropped <- sum(!macf_idx)
if (n_dropped > 0) cat("  Dropped", n_dropped, "obs with missing W_f1\n")

df_macf  <- df_f1[macf_idx, ]
X_macf   <- X_f1[macf_idx, ]
W_macf   <- droplevels(as.factor(df_macf$W_f1))
wt_macf  <- df_macf$f3bypnlwt
W_macf_num <- as.numeric(as.character(W_macf))

cat("N:", nrow(df_macf),
    "| None (0):", sum(W_macf == 0),
    "| Low (1):",  sum(W_macf == 1),
    "| High (2):", sum(W_macf == 2), "\n\n")

# ── Propensity forest (once, shared across outcomes) ─────────────────────────
cat("Fitting MACF propensity forest (W_f1 intensity, shared)...\n")
forest.W_macf <- probability_forest(as.matrix(X_macf), W_macf,
                                     sample.weights = wt_macf,
                                     num.trees = 2000, honesty = TRUE)
W.hat_macf <- predict(forest.W_macf)$predictions   # N × 3: P(W=0|X), P(W=1|X), P(W=2|X)
cat("W.hat_macf col means (should sum to 1 per row):",
    round(colMeans(W.hat_macf), 3), "\n\n")

# Overlap plot
W.hat_wide_macf <- as.data.frame(W.hat_macf)
colnames(W.hat_wide_macf) <- c("None", "Low", "High")
overlap_df_macf <- tidyr::pivot_longer(W.hat_wide_macf, everything(),
                                        names_to = "Arm", values_to = "P_hat")
overlap_df_macf$Arm <- factor(overlap_df_macf$Arm,
                              levels = c("None", "Low", "High")) #keep run order, not alphabetical
p_macf_overlap <- ggplot(overlap_df_macf, aes(x = P_hat, fill = Arm, color = Arm)) +
  geom_density(alpha = 0.6) +
  labs(fill = "Intensity", color = "Intensity",
       x = "Estimated Propensity",
       title = "Overlap: Work Intensity (12th grade)") +
  expand_limits(x = c(0, 1)) + 
  theme_minimal()
ggsave(file.path(outdir_f1, "overlap_macf.tiff"),
       p_macf_overlap, width = 6, height = 4, dpi = 600)

# ---- Combined overlap figure: 10th grade, 12th grade, 12th-grade intensity ----
# Panels A and B on the top row, panel C spanning the bottom row. The three
# standalone files above are still written; this is an additional figure.
pa <- p_overlap_by    + labs(title = NULL, subtitle = "Work (10th grade)") +
  theme(legend.position = "bottom")
pb <- p_overlap_f1    + labs(title = NULL, subtitle = "Work (12th grade)") +
  theme(legend.position = "bottom")
pc <- p_macf_overlap  + labs(title = NULL, subtitle = "Work intensity (12th grade)") +
  theme(legend.position = "bottom")

combined_overlap <- pa + pb + pc +
  plot_layout(design = "AB\nCC") +
  plot_annotation(tag_levels = "A")

ggsave("output/overlap_combined.tiff", combined_overlap,
       width = 9, height = 7, dpi = 600)
cat("Saved: output/overlap_combined.tiff\n")

# IPW per observation for balance table
ipw_macf <- numeric(nrow(df_macf))
for (k in 0:2) {
  arm_k <- W_macf_num == k
  ipw_macf[arm_k] <- 1 / W.hat_macf[arm_k, k + 1]
}
w_comb_macf <- wt_macf * ipw_macf

# ---- MACF covariate balance — pairwise (arm 1 vs 0, arm 2 vs 0) ----
X_cat_macf <- X_macf[, dummy_cols]

make_macf_balance_row <- function(v, is_binary, k, svy_base_pair, svy_ipw_pair) {
  fml    <- as.formula(paste0("~`", v, "`"))
  se_col <- "se"

  base_by <- svyby(fml, ~W_pair, svy_base_pair, svymean, na.rm = TRUE)
  ipw_by  <- svyby(fml, ~W_pair, svy_ipw_pair,  svymean, na.rm = TRUE)

  w_base <- as.numeric(as.character(base_by$W_pair))
  w_ipw  <- as.numeric(as.character(ipw_by$W_pair))

  m0_base  <- base_by[w_base == 0, v];  se0_base <- base_by[w_base == 0, se_col]
  m1_base  <- base_by[w_base == 1, v];  se1_base <- base_by[w_base == 1, se_col]
  m0_ipw   <- ipw_by[w_ipw  == 0, v];  se0_ipw  <- ipw_by[w_ipw  == 0, se_col]
  m1_ipw   <- ipw_by[w_ipw  == 1, v];  se1_ipw  <- ipw_by[w_ipw  == 1, se_col]

  overall_svy <- svymean(fml, svy_base_pair, na.rm = TRUE)
  m_overall   <- as.numeric(overall_svy)
  se_overall  <- as.numeric(SE(overall_svy))

  if (is_binary) {
    denom <- sqrt(max(m_overall * (1 - m_overall), 1e-10))
  } else {
    var_by <- svyby(fml, ~W_pair, svy_base_pair, svyvar, na.rm = TRUE)
    w_var  <- as.numeric(as.character(var_by$W_pair))
    denom  <- sqrt(max((var_by[w_var == 0, v] + var_by[w_var == 1, v]) / 2, 1e-10))
  }

  data.frame(
    comparison      = paste0("Arm ", k, " vs 0"),
    variable        = v,
    label           = get_label(v),
    mean_overall    = round(m_overall,  4),
    se_overall      = round(se_overall, 4),
    mean_treat_base = round(m1_base,    4),
    se_treat_base   = round(se1_base,   4),
    mean_ctrl_base  = round(m0_base,    4),
    se_ctrl_base    = round(se0_base,   4),
    smd_base        = round((m1_base - m0_base) / denom, 4),
    mean_treat_ipw  = round(m1_ipw,     4),
    se_treat_ipw    = round(se1_ipw,    4),
    mean_ctrl_ipw   = round(m0_ipw,     4),
    se_ctrl_ipw     = round(se0_ipw,    4),
    smd_ipw         = round((m1_ipw - m0_ipw) / denom, 4),
    row.names       = NULL
  )
}

# Build balance table for arm 1 vs 0
arm_idx_1    <- W_macf_num %in% c(0, 1)
W_pair_1     <- as.numeric(W_macf_num[arm_idx_1] == 1)
ipw_pair_1   <- ifelse(W_pair_1 == 1, 1 / W.hat_macf[arm_idx_1, 2],
                                       1 / W.hat_macf[arm_idx_1, 1])
w_pair_comb_1 <- wt_macf[arm_idx_1] * ipw_pair_1

svy_base_pair_1 <- svydesign(
  id = ~psu, strata = ~strat_id, nest = TRUE, weights = ~wt_base,
  data = cbind(df_macf[arm_idx_1, c("psu", "strat_id")],
               data.frame(W_pair = W_pair_1, wt_base = wt_macf[arm_idx_1]),
               X_cat_macf[arm_idx_1, ], df_macf[arm_idx_1, vl_cont]))
svy_ipw_pair_1 <- svydesign(
  id = ~psu, strata = ~strat_id, nest = TRUE, weights = ~w_pair_comb,
  data = cbind(df_macf[arm_idx_1, c("psu", "strat_id")],
               data.frame(W_pair = W_pair_1, w_pair_comb = w_pair_comb_1),
               X_cat_macf[arm_idx_1, ], df_macf[arm_idx_1, vl_cont]))

bal_bin_1  <- Filter(Negate(is.null), lapply(names(X_cat_macf), function(v) {
  tryCatch(make_macf_balance_row(v, is_binary = TRUE,  k = 1, svy_base_pair_1, svy_ipw_pair_1), error = function(e) NULL)
}))
bal_cont_1 <- Filter(Negate(is.null), lapply(vl_cont, function(v) {
  tryCatch(make_macf_balance_row(v, is_binary = FALSE, k = 1, svy_base_pair_1, svy_ipw_pair_1), error = function(e) NULL)
}))

# Build balance table for arm 2 vs 0
arm_idx_2    <- W_macf_num %in% c(0, 2)
W_pair_2     <- as.numeric(W_macf_num[arm_idx_2] == 2)
ipw_pair_2   <- ifelse(W_pair_2 == 1, 1 / W.hat_macf[arm_idx_2, 3],
                                       1 / W.hat_macf[arm_idx_2, 1])
w_pair_comb_2 <- wt_macf[arm_idx_2] * ipw_pair_2

svy_base_pair_2 <- svydesign(
  id = ~psu, strata = ~strat_id, nest = TRUE, weights = ~wt_base,
  data = cbind(df_macf[arm_idx_2, c("psu", "strat_id")],
               data.frame(W_pair = W_pair_2, wt_base = wt_macf[arm_idx_2]),
               X_cat_macf[arm_idx_2, ], df_macf[arm_idx_2, vl_cont]))
svy_ipw_pair_2 <- svydesign(
  id = ~psu, strata = ~strat_id, nest = TRUE, weights = ~w_pair_comb,
  data = cbind(df_macf[arm_idx_2, c("psu", "strat_id")],
               data.frame(W_pair = W_pair_2, w_pair_comb = w_pair_comb_2),
               X_cat_macf[arm_idx_2, ], df_macf[arm_idx_2, vl_cont]))

bal_bin_2  <- Filter(Negate(is.null), lapply(names(X_cat_macf), function(v) {
  tryCatch(make_macf_balance_row(v, is_binary = TRUE,  k = 2, svy_base_pair_2, svy_ipw_pair_2), error = function(e) NULL)
}))
bal_cont_2 <- Filter(Negate(is.null), lapply(vl_cont, function(v) {
  tryCatch(make_macf_balance_row(v, is_binary = FALSE, k = 2, svy_base_pair_2, svy_ipw_pair_2), error = function(e) NULL)
}))

balance_macf_low  <- do.call(rbind, c(bal_bin_1, bal_cont_1))
balance_macf_high <- do.call(rbind, c(bal_bin_2, bal_cont_2))
write.csv(balance_macf_low,  file.path(outdir_f1, "balance_macf_low.csv"),  row.names = FALSE)
write.csv(balance_macf_high, file.path(outdir_f1, "balance_macf_high.csv"), row.names = FALSE)
cat("MACF balance tables saved:", file.path(outdir_f1, "balance_macf_low.csv"),
    "and", file.path(outdir_f1, "balance_macf_high.csv"), "\n\n")

# ── Outcome: BA attainment (MACF) ────────────────────────────────────────────
cat("--- MACF Outcome: BA or Higher ---\n")
Y_macf_bach <- as.numeric(df_macf$Y_bachelors_plus)

cat("  Fitting outcome forest (MACF, bachelors)...\n")
forest.Y_macf_bach <- regression_forest(as.matrix(X_macf), Y_macf_bach,
                                         sample.weights = wt_macf,
                                         num.trees = 1500, tune.parameters = "all")
Y.hat_macf_bach <- predict(forest.Y_macf_bach)$predictions

cat("  Fitting initial MACF for variable selection (bachelors)...\n")
forest.init_macf_bach <- multi_arm_causal_forest(as.matrix(X_macf), Y_macf_bach, W_macf,
                                                   sample.weights = wt_macf,
                                                   W.hat = W.hat_macf,
                                                   Y.hat = Y.hat_macf_bach,
                                                   num.trees = 2000)
varimp_macf_bach <- variable_importance(forest.init_macf_bach)
n_keep_macf_bach <- max(1L, round(ncol(X_macf) / 5))
keep_macf_bach   <- colnames(X_macf)[order(varimp_macf_bach, decreasing = TRUE)[seq_len(n_keep_macf_bach)]]
cat("  Top fifth vars (", n_keep_macf_bach, "):", paste(keep_macf_bach, collapse = ", "), "\n")
X_sel_macf_bach <- X_macf[, keep_macf_bach]

# Base MACF (untrimmed)
res_base_bach <- fit_macf_and_save(
  Y_macf_bach, W_macf, wt_macf, W.hat_macf, Y.hat_macf_bach, X_sel_macf_bach,
  outdir = outdir_f1, out_var = "Y_bachelors_plus", out_label = "BA or Higher",
  trt_label = "Work Intensity (12th grade)", file_suffix = "", plot_suffix = "(untrimmed)")

# AIPW scores — used for both het plots and BLP
DR_scores_bach <- drop(get_scores(res_base_bach$macf))   # n x 2 after drop()

plot_macf_het(res_base_bach$macf, DR_scores_bach, wt_macf, X_sel_macf_bach,
              outdir = outdir_f1, out_var = "Y_bachelors_plus",
              out_label_str = "BA or Higher", trt_label = "Work Intensity (12th grade)",
              file_suffix = "")

# BLP via AIPW scores — one variable at a time (HC3, weighted)
varimp_macf_bach_final <- variable_importance(res_base_bach$macf)
n_blp_macf_bach    <- min(5L, ncol(X_sel_macf_bach))
blp_vars_macf_bach <- colnames(X_sel_macf_bach)[order(varimp_macf_bach_final,
                                                        decreasing = TRUE)[seq_len(n_blp_macf_bach)]]
cat("  BLP (AIPW scores, HC3, weighted, univariate) — top vars:",
    paste(blp_vars_macf_bach, collapse = ", "), "\n")

blp_macf_bach_low <- lapply(blp_vars_macf_bach, function(var) {
  A_v <- data.frame(x = X_sel_macf_bach[[var]])
  ok  <- complete.cases(cbind(DR_scores_bach[, 1], A_v))
  fit <- lm(DR_scores_bach[ok, 1] ~ x, data = A_v[ok, , drop = FALSE], weights = wt_macf[ok])
  coeftest(fit, vcov = vcovHC(fit, type = "HC3"))
})
names(blp_macf_bach_low) <- blp_vars_macf_bach

blp_macf_bach_high <- lapply(blp_vars_macf_bach, function(var) {
  A_v <- data.frame(x = X_sel_macf_bach[[var]])
  ok  <- complete.cases(cbind(DR_scores_bach[, 2], A_v))
  fit <- lm(DR_scores_bach[ok, 2] ~ x, data = A_v[ok, , drop = FALSE], weights = wt_macf[ok])
  coeftest(fit, vcov = vcovHC(fit, type = "HC3"))
})
names(blp_macf_bach_high) <- blp_vars_macf_bach

cat("  BLP: Low vs None\n")
for (var in blp_vars_macf_bach) { cat("   ", var, ":\n"); print(blp_macf_bach_low[[var]]) }
cat("  BLP: High vs None\n")
for (var in blp_vars_macf_bach) { cat("   ", var, ":\n"); print(blp_macf_bach_high[[var]]) }

blp_macf_bach_low_df <- do.call(rbind, lapply(blp_vars_macf_bach, function(var) {
  ct <- as.data.frame(unclass(blp_macf_bach_low[[var]]))
  data.frame(source_variable = var, term = rownames(ct),
             Estimate  = round(ct$Estimate,        6),
             Std_Error = round(ct[["Std. Error"]], 6),
             t_value   = round(ct[["t value"]],    4),
             p_value   = ct[["Pr(>|t|)"]],
             stringsAsFactors = FALSE, row.names = NULL)
}))
slope_bl <- blp_macf_bach_low_df$term != "(Intercept)"
blp_macf_bach_low_df$p_adj_BH <- NA_real_
blp_macf_bach_low_df$p_adj_BH[slope_bl] <- p.adjust(blp_macf_bach_low_df$p_value[slope_bl], method = "BH")
write.csv(blp_macf_bach_low_df, file = file.path(outdir_f1, "blp_macf_BA_low.csv"), row.names = FALSE)

blp_macf_bach_high_df <- do.call(rbind, lapply(blp_vars_macf_bach, function(var) {
  ct <- as.data.frame(unclass(blp_macf_bach_high[[var]]))
  data.frame(source_variable = var, term = rownames(ct),
             Estimate  = round(ct$Estimate,        6),
             Std_Error = round(ct[["Std. Error"]], 6),
             t_value   = round(ct[["t value"]],    4),
             p_value   = ct[["Pr(>|t|)"]],
             stringsAsFactors = FALSE, row.names = NULL)
}))
slope_bh <- blp_macf_bach_high_df$term != "(Intercept)"
blp_macf_bach_high_df$p_adj_BH <- NA_real_
blp_macf_bach_high_df$p_adj_BH[slope_bh] <- p.adjust(blp_macf_bach_high_df$p_value[slope_bh], method = "BH")
write.csv(blp_macf_bach_high_df, file = file.path(outdir_f1, "blp_macf_BA_high.csv"), row.names = FALSE)

cat("  MACF BLP saved:", file.path(outdir_f1, "blp_macf_BA_low.csv"),
    "and", file.path(outdir_f1, "blp_macf_BA_high.csv"), "\n")

# Trimmed MACF (sensitivity)
alpha_lo <- 0.1;  alpha_hi <- 0.9
trim_idx_bach <- apply(W.hat_macf, 1, min) >= alpha_lo & apply(W.hat_macf, 1, max) <= alpha_hi
cat("  Sensitivity trim: dropping", sum(!trim_idx_bach),
    "obs with any arm propensity outside [", alpha_lo, ",", alpha_hi, "]\n")
res_trim_bach <- fit_macf_and_save(
  Y_macf_bach[trim_idx_bach], droplevels(W_macf[trim_idx_bach]), wt_macf[trim_idx_bach],
  W.hat_macf[trim_idx_bach, ], Y.hat_macf_bach[trim_idx_bach],
  X_sel_macf_bach[trim_idx_bach, ],
  outdir = outdir_f1, out_var = "Y_bachelors_plus", out_label = "BA or Higher",
  trt_label = "Work Intensity (12th grade)", file_suffix = "_trimmed", plot_suffix = "(trimmed)")
cat("\n")

# ── Outcome: Annual earnings (MACF) ──────────────────────────────────────────
cat("--- MACF Outcome: Annual Earnings ---\n")
Y_macf_inc <- as.numeric(df_macf$Y_income)

cat("  Fitting outcome forest (MACF, income)...\n")
forest.Y_macf_inc <- regression_forest(as.matrix(X_macf), Y_macf_inc,
                                        sample.weights = wt_macf,
                                        num.trees = 1500, tune.parameters = "all")
Y.hat_macf_inc <- predict(forest.Y_macf_inc)$predictions

cat("  Fitting initial MACF for variable selection (income)...\n")
forest.init_macf_inc <- multi_arm_causal_forest(as.matrix(X_macf), Y_macf_inc, W_macf,
                                                  sample.weights = wt_macf,
                                                  W.hat = W.hat_macf,
                                                  Y.hat = Y.hat_macf_inc,
                                                  num.trees = 2000)
varimp_macf_inc <- variable_importance(forest.init_macf_inc)
n_keep_macf_inc <- max(1L, round(ncol(X_macf) / 5))
keep_macf_inc   <- colnames(X_macf)[order(varimp_macf_inc, decreasing = TRUE)[seq_len(n_keep_macf_inc)]]
cat("  Top fifth vars (", n_keep_macf_inc, "):", paste(keep_macf_inc, collapse = ", "), "\n")
X_sel_macf_inc <- X_macf[, keep_macf_inc]

# Base MACF (untrimmed)
res_base_inc <- fit_macf_and_save(
  Y_macf_inc, W_macf, wt_macf, W.hat_macf, Y.hat_macf_inc, X_sel_macf_inc,
  outdir = outdir_f1, out_var = "Y_income", out_label = "Annual Earnings ($)",
  trt_label = "Work Intensity (12th grade)", file_suffix = "", plot_suffix = "(untrimmed)")

# AIPW scores — used for both het plots and BLP
DR_scores_inc <- drop(get_scores(res_base_inc$macf))   # n x 2 after drop()

plot_macf_het(res_base_inc$macf, DR_scores_inc, wt_macf, X_sel_macf_inc,
              outdir = outdir_f1, out_var = "Y_income",
              out_label_str = "Annual Earnings ($)", trt_label = "Work Intensity (12th grade)",
              file_suffix = "")

# BLP via AIPW scores — one variable at a time (HC3, weighted)
varimp_macf_inc_final <- variable_importance(res_base_inc$macf)
n_blp_macf_inc    <- min(5L, ncol(X_sel_macf_inc))
blp_vars_macf_inc <- colnames(X_sel_macf_inc)[order(varimp_macf_inc_final,
                                                      decreasing = TRUE)[seq_len(n_blp_macf_inc)]]
cat("  BLP (AIPW scores, HC3, weighted, univariate) — top vars:",
    paste(blp_vars_macf_inc, collapse = ", "), "\n")

blp_macf_inc_low <- lapply(blp_vars_macf_inc, function(var) {
  A_v <- data.frame(x = X_sel_macf_inc[[var]])
  ok  <- complete.cases(cbind(DR_scores_inc[, 1], A_v))
  fit <- lm(DR_scores_inc[ok, 1] ~ x, data = A_v[ok, , drop = FALSE], weights = wt_macf[ok])
  coeftest(fit, vcov = vcovHC(fit, type = "HC3"))
})
names(blp_macf_inc_low) <- blp_vars_macf_inc

blp_macf_inc_high <- lapply(blp_vars_macf_inc, function(var) {
  A_v <- data.frame(x = X_sel_macf_inc[[var]])
  ok  <- complete.cases(cbind(DR_scores_inc[, 2], A_v))
  fit <- lm(DR_scores_inc[ok, 2] ~ x, data = A_v[ok, , drop = FALSE], weights = wt_macf[ok])
  coeftest(fit, vcov = vcovHC(fit, type = "HC3"))
})
names(blp_macf_inc_high) <- blp_vars_macf_inc

cat("  BLP: Low vs None\n")
for (var in blp_vars_macf_inc) { cat("   ", var, ":\n"); print(blp_macf_inc_low[[var]]) }
cat("  BLP: High vs None\n")
for (var in blp_vars_macf_inc) { cat("   ", var, ":\n"); print(blp_macf_inc_high[[var]]) }

blp_macf_inc_low_df <- do.call(rbind, lapply(blp_vars_macf_inc, function(var) {
  ct <- as.data.frame(unclass(blp_macf_inc_low[[var]]))
  data.frame(source_variable = var, term = rownames(ct),
             Estimate  = round(ct$Estimate,        6),
             Std_Error = round(ct[["Std. Error"]], 6),
             t_value   = round(ct[["t value"]],    4),
             p_value   = ct[["Pr(>|t|)"]],
             stringsAsFactors = FALSE, row.names = NULL)
}))
slope_il <- blp_macf_inc_low_df$term != "(Intercept)"
blp_macf_inc_low_df$p_adj_BH <- NA_real_
blp_macf_inc_low_df$p_adj_BH[slope_il] <- p.adjust(blp_macf_inc_low_df$p_value[slope_il], method = "BH")
write.csv(blp_macf_inc_low_df, file = file.path(outdir_f1, "blp_macf_INC_low.csv"), row.names = FALSE)

blp_macf_inc_high_df <- do.call(rbind, lapply(blp_vars_macf_inc, function(var) {
  ct <- as.data.frame(unclass(blp_macf_inc_high[[var]]))
  data.frame(source_variable = var, term = rownames(ct),
             Estimate  = round(ct$Estimate,        6),
             Std_Error = round(ct[["Std. Error"]], 6),
             t_value   = round(ct[["t value"]],    4),
             p_value   = ct[["Pr(>|t|)"]],
             stringsAsFactors = FALSE, row.names = NULL)
}))
slope_ih <- blp_macf_inc_high_df$term != "(Intercept)"
blp_macf_inc_high_df$p_adj_BH <- NA_real_
blp_macf_inc_high_df$p_adj_BH[slope_ih] <- p.adjust(blp_macf_inc_high_df$p_value[slope_ih], method = "BH")
write.csv(blp_macf_inc_high_df, file = file.path(outdir_f1, "blp_macf_INC_high.csv"), row.names = FALSE)

cat("  MACF BLP saved:", file.path(outdir_f1, "blp_macf_INC_low.csv"),
    "and", file.path(outdir_f1, "blp_macf_INC_high.csv"), "\n")

# Trimmed MACF (sensitivity) — trim_idx uses same W.hat_macf so trim set is identical
trim_idx_inc <- trim_idx_bach   # same propensities, same trim
cat("  Sensitivity trim: dropping", sum(!trim_idx_inc),
    "obs with any arm propensity outside [", alpha_lo, ",", alpha_hi, "]\n")
res_trim_inc <- fit_macf_and_save(
  Y_macf_inc[trim_idx_inc], droplevels(W_macf[trim_idx_inc]), wt_macf[trim_idx_inc],
  W.hat_macf[trim_idx_inc, ], Y.hat_macf_inc[trim_idx_inc],
  X_sel_macf_inc[trim_idx_inc, ],
  outdir = outdir_f1, out_var = "Y_income", out_label = "Annual Earnings ($)",
  trt_label = "Work Intensity (12th grade)", file_suffix = "_trimmed", plot_suffix = "(trimmed)")
cat("\n")

# ── MACF results table ────────────────────────────────────────────────────────

make_macf_ate_rows <- function(label, ate_base, ate_trim) {
  rbind(
    data.frame(Outcome    = label,
               Comparison = "Low (1-20 hrs) vs None",
               ATE_base   = round(ate_base[1, "estimate"], 4),
               SE_base    = round(ate_base[1, "std.err"],  4),
               ATE_trim   = round(ate_trim[1, "estimate"], 4),
               SE_trim    = round(ate_trim[1, "std.err"],  4),
               row.names  = NULL),
    data.frame(Outcome    = label,
               Comparison = "High (21+ hrs) vs None",
               ATE_base   = round(ate_base[2, "estimate"], 4),
               SE_base    = round(ate_base[2, "std.err"],  4),
               ATE_trim   = round(ate_trim[2, "estimate"], 4),
               SE_trim    = round(ate_trim[2, "std.err"],  4),
               row.names  = NULL)
  )
}

macf_ate_rows <- rbind(
  make_macf_ate_rows("BA or Higher",       res_base_bach$ate, res_trim_bach$ate),
  make_macf_ate_rows("Annual Earnings ($)", res_base_inc$ate,  res_trim_inc$ate)
)

write.csv(macf_ate_rows, file = file.path(outdir_f1, "macf_ate.csv"), row.names = FALSE)
cat("MACF ATE table saved:", file.path(outdir_f1, "macf_ate.csv"), "\n")
cat("MACF BLP saved to blp_macf_BA_low/high.csv and blp_macf_INC_low/high.csv\n\n")

# ==============================================================================
# Done
# ==============================================================================

cat("CF_analysis.R completed:", format(Sys.time()), "\n")
sink()
sink(type = "message")
close(msg_con)
