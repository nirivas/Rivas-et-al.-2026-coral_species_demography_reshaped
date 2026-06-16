

# Load Libraries ----------------------------------------------------------

library(tidyverse)
library(ipmr)
library(performance)
library(truncnorm)
library(MASS)
library(gridExtra)
library(ggeffects)
library(ggpubr)
library(stringr)
library(ggrepel)
library(glmmTMB)
library(png)
library(grid)
library(parallel)
library(pbapply) 



# Load Data ---------------------------------------------------------------


df_diploria = read_csv("data/df_diploria.csv") |> 
  mutate(transition = factor(transition, levels = c("t1", "t2")))


# Dlab  -------------------------------------------------------------------


## Vital Rates -------------------------------------------------------------


### Growth ------------------------------------------------------------------



g_global2_dlab = lm(log_size_next ~ log_size * transition, data = df_diploria)

summary(g_global2_dlab) 

plot(ggpredict(g_global2_dlab, terms = c("log_size", "transition")))

grow_sd_dlab_global  <- sd(resid(g_global2_dlab))


pred_dlab <- ggpredict(g_global2_dlab, terms = c("log_size", "transition"))
pred_df_dlab <- as.data.frame(pred_dlab)
# columns: x, predicted, conf.low, conf.high, group (transition)

# Relabel transitions if you want nicer labels
pred_df_dlab$group <- dplyr::recode(pred_df_dlab$group,
                                    "t1" = "2021–2022",
                                    "t2" = "2022–2023")

growth_dlab = ggplot(pred_df_dlab,
                     aes(x = x,
                         y = predicted,
                         colour = group,
                         fill   = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.5,
              colour = NA) +
  geom_line(linewidth = 1) +
  scale_fill_viridis_d(option = "H", begin = .2, end = 0.9) +
  scale_colour_viridis_d(option = "H", begin = .2, end = 0.9) +
  theme_classic() +
  labs(x     = "Size at time t",
       y     = "Size at time t + 1",
       colour = "Transition",
       fill   = "Transition") +
  theme_classic()+
  theme(axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text  = element_text(size = 14),
        legend.text = element_text(size = 12))


growth_dlab



### Survival ----------------------------------------------------------------



s_global_dlab  <- glm(survival ~ log_size * transition, family = binomial, data = df_diploria) 


summary(s_global_dlab)

plot(ggpredict(s_global_dlab, terms = c("log_size[all]", "transition")))

pred_survival_dlab <- ggpredict(s_global_dlab, terms = c("log_size[all]", "transition"))
pred_survival_df_dlab <- as.data.frame(pred_survival_dlab)

pred_survival_df_dlab$group <- dplyr::recode(pred_survival_df_dlab$group,
                                             "t1" = "2021–2022",
                                             "t2" = "2022–2023")

survival_dlab = ggplot(pred_survival_df_dlab,
                       aes(x = x,
                           y = predicted*100,
                           colour = group,
                           fill   = group)) +
  geom_ribbon(aes(ymin = conf.low*100, ymax = conf.high*100),
              alpha = 0.5,
              colour = NA) +
  geom_line(linewidth = 1) +
  scale_fill_viridis_d(option = "H", begin = .2, end = 0.9) +
  scale_colour_viridis_d(option = "H", begin = .2, end = 0.9) +
  scale_y_continuous(
    name   = "Predicted survival (%)",
    limits = c(30, 100)) +
  labs(x     = "Size at time t",
       colour = "Transition",
       fill   = "Transition") +
  theme_classic()+
  theme(axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text  = element_text(size = 14),
        legend.text = element_text(size = 12))

survival_dlab


## IPM ---------------------------------------------------------------------


### Parameters --------------------------------------------------------------



params_global1_dlab <- list(
  grow_mod   = use_vr_model(g_global2_dlab),   # lm with log_size * transition
  surv_mod   = s_global_dlab,    # glm with log_size * transition
  grow_sd    = grow_sd_dlab_global,
  transition = "t1"
) 

params_global2_dlab <- list(
  grow_mod   = use_vr_model(g_global2_dlab),   # lm with log_size * transition
  surv_mod   = s_global_dlab,      # glm with log_size * transition
  grow_sd    = grow_sd_dlab_global,
  transition = "t2"
)



### IPM Function ------------------------------------------------------------



U_dlab <- max(c(df_diploria$log_size, df_diploria$log_size_next), na.rm = TRUE) * 1.2 


L_dlab <- min(c(df_diploria$log_size, df_diploria$log_size_next), na.rm = TRUE) * 0.8



#### Transition 1 ------------------------------------------------------------


ipm_global1_dlab <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "det") %>%
  
  # --- P Kernel (Survival * Growth) ---
  define_kernel(
    name    = "P",
    family  = "CC",
    formula = s * g,
    
    s = predict(surv_mod,
                newdata = data.frame(
                  log_size   = sa_1 ,
                  transition = rep(transition, length(sa_1))
                ),
                type = "response"),
    
    g_mu = predict(grow_mod,
                   newdata = data.frame(
                     log_size   = sa_1,
                     transition = rep(transition, length(sa_1))
                   ),
                   type = "response"),
    
    g = dnorm(sa_2,g_mu, grow_sd),
    
    states        = list(c("sa")),
    data_list     = params_global1_dlab,
    uses_par_sets  = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions(fun = "norm",
                                            target = "g")
  )  %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P"),
      int_rule     = rep('midpoint'),
      state_start    = rep("sa"),
      state_end      = rep("sa")
    )
  ) %>%
  define_domains(
    sa = c(L_dlab,
           U_dlab,
           200)
  ) %>%
  define_pop_state(n_sa = runif(200)) %>%
  make_ipm(iterate = TRUE,
           iterations = 1000)



#### Transition 2 ------------------------------------------------------------



ipm_global2_dlab <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "det") %>%
  
  # --- P Kernel (Survival * Growth) ---
  define_kernel(
    name    = "P",
    family  = "CC",
    formula = s * g,
    
    s = predict(surv_mod,
                newdata = data.frame(
                  log_size   = sa_1 ,
                  transition = rep(transition, length(sa_1))
                ),
                type = "response"),
    
    g_mu = predict(grow_mod,
                   newdata = data.frame(
                     log_size   = sa_1,
                     transition = rep(transition, length(sa_1))
                   ),
                   type = "response"),
    
    g = dnorm(sa_2,g_mu, grow_sd),
    
    states        = list(c("sa")),
    data_list     = params_global2_dlab,
    uses_par_sets  = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions(fun = "norm",
                                            target = "g")
  )  %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P"),
      int_rule     = rep('midpoint'),
      state_start    = rep("sa"),
      state_end      = rep("sa")
    )
  ) %>%
  define_domains(
    sa = c(L_dlab,
           U_dlab,
           200)
  ) %>%
  define_pop_state(n_sa = runif(200)) %>%
  make_ipm(iterate = TRUE,
           iterations = 1000)




#### Lambda ------------------------------------------------------------------


lambda(ipm_global1_dlab)
ipmr::is_conv_to_asymptotic(ipm_global1_dlab)

lambda(ipm_global2_dlab)




### Deterministic Projection ------------------------------------------------



tmax <- 50

# lambdas
lam1_dlab <- lambda(ipm_global1_dlab)
lam2_dlab <- lambda(ipm_global2_dlab) 

# starting population size
N0 <- 100  # or 1 if you want relative trajectories

# make projections
proj1_dlab <- tibble(
  time = 0:tmax,
  pop_size = N0 * lam1_dlab^(0:tmax)
)

proj2_dlab <- tibble(
  time = 0:tmax,
  pop_size = N0 * lam2_dlab^(0:tmax)
)

# plot model 1
a_dlab = ggplot(proj1_dlab, aes(x = time, y = pop_size)) +
  geom_line(color = "steelblue", linewidth = 1.2) +
  labs(
    title = paste0("Pre die-off (λ = ", round(lam1_dlab, 3), ")"),
    x = "Years",
    y = "Population size"
  ) +
  theme_classic() +
  theme(axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text  = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(face = "bold"))

# plot model 2
b_dlab = ggplot(proj2_dlab, aes(x = time, y = pop_size)) +
  geom_line(color = "forestgreen", linewidth = 1.2) +
  labs(
    title = paste0("Post die-off (λ = ", round(lam2_dlab, 3), ")"),
    x = "Years",
    y = "Population size"
  ) +
  theme_classic() +
  theme(axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text  = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(face = "bold"))



# combine plots

ggarrange(a_dlab, b_dlab,
          ncol =2, nrow = 1,
          labels = c("A", "B"))



## Lambda Bootstrap --------------------------------------------------------




### Function ###


# build_ipm_lambda_dlab <- function(df_boot, transition = "t1", mesh_size = 200) {
#   
#   g_mod_boot   <- update(g_global2_dlab, data = df_boot)
#   s_mod_boot   <- update(s_global_dlab,  data = df_boot)
#   grow_sd_boot <- sd(resid(g_mod_boot), na.rm = TRUE)
#   
#   params_boot <- list(
#     grow_mod   = use_vr_model(g_mod_boot),
#     surv_mod   = s_mod_boot,
#     grow_sd    = grow_sd_boot,
#     transition = transition
#   )
#   
#   L_boot <- min(c(df_diploria$log_size, df_diploria$log_size_next), na.rm = TRUE) * 0.8
#   U_boot <- max(c(df_diploria$log_size, df_diploria$log_size_next), na.rm = TRUE) * 1.2
#   
#   ipm_boot <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "det") %>%
#     define_kernel(
#       name = "P", family = "CC", formula = s * g,
#       s    = predict(surv_mod,
#                      newdata = data.frame(log_size  = sa_1,
#                                           transition = rep(transition, length(sa_1))),
#                      type = "response"),
#       g_mu = predict(grow_mod,
#                      newdata = data.frame(log_size  = sa_1,
#                                           transition = rep(transition, length(sa_1))),
#                      type = "response"),
#       g    = dnorm(sa_2, g_mu, grow_sd),
#       states        = list(c("sa")),
#       data_list     = params_boot,
#       uses_par_sets = FALSE,
#       evict_cor     = TRUE,
#       evict_fun     = truncated_distributions(fun = "norm", target = "g")
#     ) %>%
#     define_impl(make_impl_args_list(
#       kernel_names = "P",
#       int_rule     = "midpoint",
#       state_start  = "sa",
#       state_end    = "sa"
#     )) %>%
#     define_domains(sa = c(L_boot, U_boot, mesh_size)) %>%
#     define_pop_state(n_sa = runif(mesh_size)) %>%
#     make_ipm(iterate = TRUE, iterations = 1000)
#   
#   lambda(ipm_boot)
# }
# 
# # ── One bootstrap replicate returning both t1 and t2 ──────────────────────────
# one_boot_dlab <- function(i) {
#   set.seed(i)
#   idx     <- sample(seq_len(nrow(df_diploria)), replace = TRUE)
#   df_boot <- df_diploria[idx, ]
#   
#   lam_t1 <- tryCatch(suppressWarnings(build_ipm_lambda_dlab(df_boot, transition = "t1")),
#                      error = function(e) NA_real_)
#   lam_t2 <- tryCatch(suppressWarnings(build_ipm_lambda_dlab(df_boot, transition = "t2")),
#                      error = function(e) NA_real_)
#   
#   c(t1 = lam_t1, t2 = lam_t2)
# }
# 
# # ── Windows parallel setup ─────────────────────────────────────────────────────
# n_cores <- detectCores() - 6
# cl      <- makeCluster(n_cores)
# 
# clusterExport(cl, varlist = c("df_diploria",
#                               "build_ipm_lambda_dlab",
#                               "g_global2_dlab",
#                               "s_global_dlab",
#                               "use_vr_model"))
# 
# clusterEvalQ(cl, {
#   library(ipmr)
#   library(dplyr)
# })
# 
# # ── Run with progress bar ──────────────────────────────────────────────────────
# set.seed(123)
# cat("Starting bootstrap — running 1200 replicates across", n_cores, "cores\n")
# t_start  <- proc.time()
# 
# results  <- pblapply(1:1200, one_boot_dlab, cl = cl)   # <-- progress bar here
# 
# t_end    <- proc.time()
# stopCluster(cl)
# cat("Done! Elapsed time:", round((t_end - t_start)["elapsed"] / 60, 1), "minutes\n")
# 
# # ── Parse results ──────────────────────────────────────────────────────────────
# results_df     <- do.call(rbind, results)
# 
# lambda_boot_t1 <- na.omit(results_df[, "t1.lambda"])[1:1000]
# lambda_boot_t2 <- na.omit(results_df[, "t2.lambda"])[1:1000]
# 
# cat("Valid t1 replicates:", sum(is.finite(results_df[, "t1.lambda"])), "\n")
# cat("Valid t2 replicates:", sum(is.finite(results_df[, "t2.lambda"])), "\n")
# 
# quantile(lambda_boot_t1, c(0.025, 0.975))
# quantile(lambda_boot_t2, c(0.025, 0.975))
# 
# # Point estimates
# mean_t1 <- mean(lambda_boot_t1)
# mean_t2 <- mean(lambda_boot_t2)
# 
# ci_t1 <- quantile(lambda_boot_t1, c(0.025, 0.5, 0.975))
# ci_t2 <- quantile(lambda_boot_t2, c(0.025, 0.5, 0.975))
# 
# mean_t1; ci_t1
# mean_t2; ci_t2
# 
# # Histograms
# par(mfrow = c(1, 2))
# hist(lambda_boot_t1, breaks = 30,
#      main = expression("Bootstrap " * lambda * " 2021-2022"),
#      xlab = expression(lambda))
# abline(v = ci_t1[c(1,3)], col = "red",  lty = 2)
# abline(v = ci_t1[2],       col = "blue", lwd = 2)
# 
# hist(lambda_boot_t2, breaks = 30,
#      main = expression("Bootstrap " * lambda * " 2022-2023"),
#      xlab = expression(lambda))
# abline(v = ci_t2[c(1,3)], col = "red",  lty = 2)
# abline(v = ci_t2[2],       col = "blue", lwd = 2)
# par(mfrow = c(1, 1))
# 
# # Combined ggplot
# all_lambdas_dlab <- data.frame(
#   lambda     = c(lambda_boot_t1, lambda_boot_t2),
#   transition = rep(c("2021–2022", "2022–2023"), each = 1000)
# )
# 
# write_csv(all_lambdas_dlab,"Markdowns/data/dlab_gloval_ipm_lambda_bootstrap.csv")


### Boostrap data -----------------------------------------------------------


all_lambdas_dlab = read_csv("data/dlab_global_ipm_lambda_bootstrap.csv")  

### Density Plot ###

boot_dlab = ggplot(all_lambdas_dlab, aes(x = lambda, fill = transition)) +
  geom_density(alpha = 0.35) +
  labs(
    x = expression(lambda),
    y = "Density", fill = "Transition") +
  scale_fill_viridis_d(option = "H", begin = .9, end = 0.2) +
  scale_x_continuous(limits = c(0.6, 1.1)) +
  theme_classic()+
  theme(axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text  = element_text(size = 14),
        legend.text = element_text(size = 12))
boot_dlab

boot <- all_lambdas_dlab %>%
  dplyr::select(transition, lambda) %>%
  filter(!is.na(lambda))


### Stochastic lambda projectoin --------------------------------------------

n_sims = 100
years <- 2025:2075
n_years = 51

sim <- boot %>%
  group_by(transition) %>%
  summarise(lambda_draws = list(lambda), .groups = "drop") %>%
  tidyr::uncount(n_sims, .id = "sim_id") %>%
  rowwise() %>%
  mutate(
    lambda_path = list(sample(lambda_draws, size = n_years - 1, replace = TRUE)),
    year = list(years),
    N = list(c(
      N0,
      N0 * exp(cumsum(log(lambda_path)))
    ))
  ) %>%
  unnest(c(year, N)) %>%
  ungroup()

# Plot
f2b = ggplot(sim, aes(x = year, y = N,
                      group = interaction(transition, sim_id),
                      color = transition)) +
  # annotation_custom(
  #   grob_imgdlab,
  #   xmin = 2067, xmax = 2072,   # <-- adjust to your x scale
  #   ymin = 130,
  #   ymax = 160
  # )+
  geom_line(alpha = .1, linewidth = 0.3) +
  labs(
    x = "Year",
    y = "Population size",
    color = "Transition"
  ) +
  stat_summary(
    aes(group = transition),
    fun = mean,
    geom = "line",
    linewidth = 1.2)+
  scale_color_viridis_d(option = "H",begin = .2, end = 0.9) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    axis.text  = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
f2b



## HeatMaps ---------------------------------------------------------


### P kernel ----------------------------------------------------------------



#### Transition 1 ------------------------------------------------------------


Kt1 <- ipm_global1_dlab$sub_kernels$P

wt1 <- Re(eigen(Kt1)$vectors[, 1])          # right eigenvector
vt1 <- Re(eigen(t(Kt1))$vectors[, 1])       # left eigenvector

# Ensure positive orientation
wt1 <- abs(wt1)
vt1 <- abs(vt1)

mesh_info_dlabt1<-int_mesh(ipm_global1_dlab) 

d_sat1 <- mesh_info_dlabt1$d_sa

sens_mat_dlabt1 <- outer(vt1, wt1) / sum(vt1 * wt1 * d_sat1)
elas_mat_dlabt1 <- sens_mat_dlabt1 * Kt1 / lambda(ipm_global1_dlab)

sens_df_dlabt1 <- ipm_to_df(sens_mat_dlabt1) 
elas_df_dlabt1 <- ipm_to_df(elas_mat_dlabt1)

def_theme <- theme( 
  panel.background = element_blank(),
  axis.text = element_text(size = 16), 
  axis.ticks = element_line(linewidth  = 1.5), 
  axis.ticks.length = unit(0.08, "in"), 
  axis.title.x = element_text( 
    size = 20, 
    margin = margin( 
      t = 10, 
      r = 0, 
      l = 0, 
      b = 2 
    ) 
  ),
  axis.title.y = element_text( 
    size = 20, 
    margin = margin( 
      t = 0, 
      r = 10, 
      l = 2, 
      b = 0 
    )
  ), 
  legend.text = element_text(size = 16)
)


p_df_dlabt1 <- data.frame(
  t     = mesh_info_dlabt1$sa_1,
  t_1   = mesh_info_dlabt1$sa_2,
  value = as.vector(ipm_global1_dlab$sub_kernels$P)
)

p_plt_dlabt1 <- ggplot(p_df_dlabt1) + 
  geom_tile(aes(x = t, 
                y = t_1, 
                fill = value)) + 
  geom_contour(aes(x = t, 
                   y = t_1, 
                   z = value),
               color = "black", 
               linewidth = 0.7, 
               bins = 5) + 
  scale_fill_gradient("Value", 
                      low = "red", 
                      high = "yellow") + 
  scale_x_continuous(name = "Area at t (cm²)",
                     breaks = log(c(5, 25, 100, 400, 800)),
                     labels = c(5, 25, 100, 400, 800))+
scale_y_continuous(name = "Area at t+1 (cm²)",
                   breaks = log(c(5, 25, 100, 400, 800)),
                   labels = c(5, 25, 100, 400, 800))+
  geom_abline(slope = 1, intercept = 0, color = "blue", linewidth = 1) +
  def_theme +
  theme_classic() +
  theme(axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text  = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(face = "bold")) + 
  ggtitle("P kernel Transition 1")

p_plt_dlabt1

# 1) Conditional *mean* size at t+1 for each t
ridge_mean_dlabt1 <- p_df_dlabt1 |>
  dplyr::group_by(t) |>
  dplyr::summarise(t1_bar = sum(t_1 * value) / sum(value), .groups = "drop") |>
  dplyr::mutate(increment = t1_bar - t)


# 2) Overlay that mean curve
p_plt_dlabt1 +
  geom_line(data = ridge_mean_dlabt1, aes(x = t, y = t1_bar),
            linewidth = 1.1, color = "white")




#### Transition 2 ------------------------------------------------------------


Kt2 <- ipm_global2_dlab$sub_kernels$P

wt2 <- Re(eigen(Kt2)$vectors[, 1])          # right eigenvector
vt2 <- Re(eigen(t(Kt2))$vectors[, 1])       # left eigenvector

# Ensure positive orientation
wt2 <- abs(wt2)
vt2 <- abs(vt2)

mesh_info_dlabt2<-int_mesh(ipm_global2_dlab) 

d_sat2 <- mesh_info_dlabt2$d_sa

sens_mat_dlabt2 <- outer(vt2, wt2) / sum(vt2 * wt2 * d_sat2)
elas_mat_dlabt2 <- sens_mat_dlabt2 * Kt2 / lambda(ipm_global2_dlab)

sens_df_dlabt2 <- ipm_to_df(sens_mat_dlabt2) 
elas_df_dlabt2 <- ipm_to_df(elas_mat_dlabt2)

def_theme <- theme( 
  panel.background = element_blank(),
  axis.text = element_text(size = 16), 
  axis.ticks = element_line(linewidth  = 1.5), 
  axis.ticks.length = unit(0.08, "in"), 
  axis.title.x = element_text( 
    size = 20, 
    margin = margin( 
      t = 10, 
      r = 0, 
      l = 0, 
      b = 2 
    ) 
  ),
  axis.title.y = element_text( 
    size = 20, 
    margin = margin( 
      t = 0, 
      r = 10, 
      l = 2, 
      b = 0 
    )
  ), 
  legend.text = element_text(size = 16)
)


p_df_dlabt2 <- data.frame(
  t     = mesh_info_dlabt2$sa_1,
  t_1   = mesh_info_dlabt2$sa_2,
  value = as.vector(ipm_global2_dlab$sub_kernels$P)
)

p_plt_dlabt2 <- ggplot(p_df_dlabt2) + 
  geom_tile(aes(x = t, 
                y = t_1, 
                fill = value)) + 
  geom_contour(aes(x = t, 
                   y = t_1, 
                   z = value),
               color = "black", 
               linewidth = 0.7, 
               bins = 5) + 
  scale_fill_gradient("Value", 
                      low = "red", 
                      high = "yellow") + 
  scale_x_continuous(name = "Area at t (cm²)",
                     breaks = log(c(5, 25, 100, 400, 800)),
                     labels = c(5, 25, 100, 400, 800))+
  scale_y_continuous(name = "Area at t+1 (cm²)",
                     breaks = log(c(5, 25, 100, 400, 800)),
                     labels = c(5, 25, 100, 400, 800))+
  geom_abline(slope = 1, intercept = 0, color = "blue", linewidth = 1) +
  def_theme +
  theme_classic() +
  theme(axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text  = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(face = "bold")) + 
  ggtitle("P kernel Transition 2")

p_plt_dlabt2

# 1) Conditional *mean* size at t+1 for each t
ridge_mean_dlabt2 <- p_df_dlabt2 |>
  dplyr::group_by(t) |>
  dplyr::summarise(t1_bar = sum(t_1 * value) / sum(value), .groups = "drop") |>
  dplyr::mutate(increment = t1_bar - t)


# 2) Overlay that mean curve
p_plt_dlabt2 +
  geom_line(data = ridge_mean_dlabt2, aes(x = t, y = t1_bar),
            linewidth = 1.1, color = "white")



### Elasticity --------------------------------------------------------------


#### Transition 1 ------------------------------------------------------------


sa_mesht1 <- mesh_info_dlabt1$sa

elas_df_dlabt1 <- data.frame(
  t   = mesh_info_dlabt1$sa_1,
  t_1 = mesh_info_dlabt1$sa_2,
  value = as.vector(elas_mat_dlabt1)
)

elas_plt_dlabt1 <- ggplot(elas_df_dlabt1) + 
  geom_tile(aes(x = t,
                y = t_1,
                fill = value)) + 
  geom_contour(aes(x = t,
                   y = t_1,
                   z = value), 
               color = "black",
               linewidth = 0.7, 
               bins = 5) + 
  scale_fill_gradient("Value",
                      low = "red", 
                      high = "yellow") +
  scale_x_continuous(name = "Area at t (cm²)",
                     breaks = log(c(5, 25, 50, 100, 200, 400, 800)),
                     labels = c(5, 25, 50, 100, 200, 400, 800),
                     limits = c(L_dlab, U_dlab))+
scale_y_continuous(name = "Area at t+1 (cm²)",
                   breaks = log(c(5, 25, 50, 100, 200, 400, 800)),
                   labels = c(5, 25, 50, 100, 200, 400, 800),
                   limits = c(L_dlab, U_dlab))+
  
  theme_classic() +
  theme(axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text  = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(face = "bold")) +
  ggtitle("K Elasticity Transition 1")


elas_plt_dlabt1




#### Transition 2 ------------------------------------------------------------




sa_mesht2 <- mesh_info_dlabt2$sa

elas_df_dlabt2 <- data.frame(
  t   = mesh_info_dlabt2$sa_1,
  t_1 = mesh_info_dlabt2$sa_2,
  value = as.vector(elas_mat_dlabt2)
)

elas_plt_dlabt2 <- ggplot(elas_df_dlabt2) + 
  geom_tile(aes(x = t,
                y = t_1,
                fill = value)) + 
  geom_contour(aes(x = t,
                   y = t_1,
                   z = value), 
               color = "black",
               linewidth = 0.7, 
               bins = 5) + 
  scale_fill_gradient("Value",
                      low = "red", 
                      high = "yellow") +
  scale_x_continuous(name = "Area at t (cm²)",
                     breaks = log(c(5, 25, 50, 100, 200, 400, 800)),
                     labels = c(5, 25, 50, 100, 200, 400, 800),
                     limits = c(L_dlab, U_dlab))+
  scale_y_continuous(name = "Area at t+1 (cm²)",
                     breaks = log(c(5, 25, 50, 100, 200, 400, 800)),
                     labels = c(5, 25, 50, 100, 200, 400, 800),
                     limits = c(L_dlab, U_dlab))+
  
  theme_classic() +
  theme(axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text  = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(face = "bold")) +
  ggtitle("K Elasticity Transition 2")


elas_plt_dlabt2



