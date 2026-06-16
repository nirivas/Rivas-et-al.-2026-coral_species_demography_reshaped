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


df_porites = read_csv("Markdowns/data/processed/df_porites.csv") |> 
  mutate(transition = factor(transition, levels = c("t1", "t2")))


# Porites -----------------------------------------------------------------



## Size Frquency Dist ------------------------------------------------------


df_long = df_porites |> 
  pivot_longer(cols = c(area1, area2),
               names_to = "area_type",
               values_to = "size")

df_longt1 = df_long |> 
  filter(transition == "t1")

df_longt2 = df_long |> 
  filter(transition == "t2")


### Plots -------------------------------------------------------------------



## Plot t1 ##

pasfd1 = ggplot(df_longt1, aes(size, fill = area_type)) +
  geom_density(alpha = 0.9) +
  scale_x_continuous(
    name   = "Colony size (cm²)",
    limits = c(0, 250),
    expand = c(.01, 0)
  ) +
  scale_y_continuous(
    name   = "Density",
    limits = c(0, 0.041),   # adjust max as needed
    expand = c(0, .0007)
  ) +
  scale_fill_viridis_d(
    option = "H",
    begin = .2,
    end = .9,
    labels = c(
      area1 = "t",
      area2 = "t + 1")) +
  theme_classic()+
  labs(fill = "Time",
       title = "**2021 - 2022**")+
  theme(axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text  = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = ggtext::element_markdown(
          hjust = 0.5,
          size = 16))

## Plot t2 ##

pasfd2 = ggplot(df_longt2, aes(size, fill = area_type)) +
  geom_density(alpha = 0.9) +
  scale_x_continuous(
    name   = "Colony size (cm²)",
    limits = c(0, 250),
    expand = c(.01, 0)
  ) +
  scale_y_continuous(
    name   = "Density",
    limits = c(0, 0.041),   # adjust max as needed
    expand = c(0, .0007)
  ) +
  scale_fill_viridis_d(
    option = "H",
    begin = .2,
    end = .9,
    labels = c(
      area1 = "t",
      area2 = "t + 1")) +
  theme_classic()+
  labs(fill = "Time",
       title = "**2022 - 2023**")+
  theme(axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text  = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = ggtext::element_markdown(
          hjust = 0.5,
          size = 16)) 


#ggsave('figures/results/sup/SFD.png', units = 'in', width = 15, height = 12, dpi = 600)

ggpubr::ggarrange(
  pasfd1, pasfd2,
  ncol = 2, nrow = 2,
  labels = c("A)", "B)", "C)", "D)"),
  common.legend = TRUE, legend = "bottom")

#ggsave('figures/results/sup/SFD_all.png', units = 'in', width = 14, height = 12, dpi = 600)


### Stats -------------------------------------------------------------------

stats_summary <- df_long %>%
  group_by(transition, area_type) %>%
  summarise(
    n = sum(!is.na(size)),
    mean_size = mean(size, na.rm = TRUE),
    median_size = median(size, na.rm = TRUE),
    sd_size = sd(size, na.rm = TRUE),
    min_size = min(size, na.rm = TRUE),
    max_size = max(size, na.rm = TRUE),
    p10 = quantile(size, 0.10, na.rm = TRUE),
    p90 = quantile(size, 0.90, na.rm = TRUE)
  )

stats_summary



ks_results <- df_long %>%
  group_by(transition) %>%
  summarise(
    ks_p = ks.test(size[area_type == "area1"],
                   size[area_type == "area2"])$p.value
  )

ks_results


## 2. Mean + median size for each year (area_type) within each transition
summary_stats <- df_long %>%
  group_by(transition, area_type) %>%
  summarise(
    mean_size   = mean(size, na.rm = TRUE),
    median_size = median(size, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = area_type,
    values_from = c(mean_size, median_size),
    names_glue  = "{.value}_{area_type}"
  )

## 3. Join stats + KS p, format, and relabel transitions
ks_table_data <- summary_stats %>%
  left_join(ks_results, by = "transition") %>%
  mutate(
    transition = recode(transition,
                        t1 = "2021–2022",
                        t2 = "2022–2023"),
    across(starts_with("mean_size"),   ~ round(., 2)),
    across(starts_with("median_size"), ~ round(., 2)),
    ks_p = pvalue(ks_p, accuracy = 0.001)
  )

## 4. Convert to flextable
ks_flex <- ks_table_data %>%
  flextable() %>%
  set_caption("Table X. Kolmogorov–Smirnov tests and summary statistics for colony size distributions between years") %>%
  set_header_labels(
    transition          = "Transition",
    mean_size_area1     = "Mean size (previous year)",
    mean_size_area2     = "Mean size (next year)",
    median_size_area1   = "Median size (previous year)",
    median_size_area2   = "Median size (next year)",
    ks_p                = "KS p-value"
  ) %>%
  autofit()

## 5. Add to Word doc and save
doc <- read_docx() %>%
  body_add_par("Kolmogorov–Smirnov tests for size distributions", 
               style = "heading 1") %>%
  body_add_flextable(ks_flex)

print(doc, target = "Markdowns/tables/processed/ks_size_distributions.docx")

df_small <- df_long %>%
  filter(action == "born",
         area_type == "area2",
         size < 5)



ggplot(df_small, aes(x = size, fill = transition)) +
  geom_histogram(binwidth = 0.5, alpha = 0.8, position = "identity") +
  scale_fill_viridis_d(option = "H", begin = .9,end = .2,
                       labels = c(t1 = "2022", t2 = "2023")) +
  scale_x_continuous(
    name   = "Colony size (mm²)",
    limits = c(0, 5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name   = "Count",
    limits = c(0, 180),   # adjust max as needed
    expand = c(0, 0)
  ) +
  labs(fill = "Recruits in")+
  theme_classic()+ 
  theme(axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text  = element_text(size = 14),
        legend.text = element_text(size = 12)) 

#ggsave('figures/results/sup/Recruits_bar.png', units = 'in', width = 10, height = 8, dpi = 600)


## 1) Restrict to P. astreoides if df has species column
# df_pa <- df %>% filter(species == "P. astreoides")

## 2) Pivot as you already do
df_long <- df_porites |> 
  pivot_longer(cols = c(area1, area2),
               names_to = "area_type",
               values_to = "size")

## 3) Compare across transitions using the SAME area_type (area2)
df_pa_area2 <- df_long %>%
  filter(area_type == "area2") %>%            # compare "next year" only
  filter(!is.na(size))

ks_t1_vs_t2_area2 <- ks.test(
  df_pa_area2$size[df_pa_area2$transition == "t1"],  # 2022 sizes
  df_pa_area2$size[df_pa_area2$transition == "t2"]   # 2023 sizes
)

ks_t1_vs_t2_area2$p.value

ggplot(df_pa_area2, aes(x = size, fill = transition)) +
  geom_density(alpha = 0.8) +
  scale_fill_viridis_d(
    option = "H",
    begin = .9, end = .2,
    labels = c(t1 = "2022 (end of 2021–2022)", t2 = "2023 (end of 2022–2023)")
  ) +
  scale_x_continuous(name = "Colony size (cm²)", limits = c(0, 250), expand = c(.01, 0)) +
  scale_y_continuous(name = "Density", expand = c(0, .0007)) +
  theme_classic() +
  labs(fill = "Transition outcome year")


## Vital Rates -------------------------------------------------------------


### Growth  ------------------------------------------------------------


g_global = lm(log_size_next ~ log_size, data = df_porites)
g_global2 = lm(log_size_next ~ log_size * transition, data = df_porites)

compare_performance(g_global, g_global2)

summary(g_global2)

plot(ggpredict(g_global2, terms = c("log_size", "transition")))

pred <- ggpredict(g_global2, terms = c("log_size", "transition"))
pred_df <- as.data.frame(pred)
# columns: x, predicted, conf.low, conf.high, group (transition)

# Relabel transitions if you want nicer labels
pred_df$group <- dplyr::recode(pred_df$group,
                               "t1" = "2021–2022",
                               "t2" = "2022–2023")

growth = ggplot(pred_df,
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

#ggsave('figures/results/Growth_model.png', units = 'in', width = 10, height = 8, dpi = 600)

# Tidy the lm model
# g_tidy <- tidy(g_global2)
# 
# # Format and turn into flextable
# g_table <- g_tidy %>%
#   mutate(
#     term = recode(term,
#                   "(Intercept)"          = "Intercept",
#                   "log_size"             = "log(Size)",
#                   "transitiont2"         = "Transition (t2 vs t1)",
#                   "log_size:transitiont2"= "log(Size) × Transition"),
#     across(c(estimate, std.error, statistic), ~ round(., 3)),
#     p.value = pvalue(p.value, accuracy = 0.001)
#   ) %>%
#   flextable() %>%
#   set_caption("Table X. Linear model of log colony size at time t+1 as a function of log size and transition") %>%
#   set_header_labels(
#     term       = "Term",
#     estimate   = "Estimate",
#     std.error  = "Std. Error",
#     statistic  = "t value",
#     p.value    = "p-value"
#   ) %>%
#   autofit()
# 
# # Create Word doc and add the table
# doc <- read_docx() %>%
#   body_add_par("Linear model of growth (g_global2)", style = "heading 1") %>%
#   body_add_flextable(g_table)
# 
# # Save the Word document
# print(doc, target = "tables/growth_global_lm_table.docx")



### Survival ----------------------------------------------------------------


s_global  <- glm(survival ~ log_size * transition, family = binomial, data = df_porites)


summary(s_global)

plot(ggpredict(s_global, terms = c("log_size[all]", "transition")))

pred_survival <- ggpredict(s_global, terms = c("log_size[all]", "transition"))
pred_survival_df <- as.data.frame(pred_survival)

pred_survival_df$group <- dplyr::recode(pred_survival_df$group,
                                        "t1" = "2021–2022",
                                        "t2" = "2022–2023")

survival = ggplot(pred_survival_df,
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
    limits = c(65, 100)) +
  labs(x     = "Size at time t",
       colour = "Transition",
       fill   = "Transition") +
  theme_classic()+
  theme(axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text  = element_text(size = 14),
        legend.text = element_text(size = 12))
# 
# ggsave('figures/results/survival_model.png', units = 'in', width = 10, height = 8, dpi = 600)

# s_tidy <- tidy(s_global)
# 
# # 2. Format table
# s_table <- s_tidy %>%
#   mutate(
#     term = recode(term,
#                   "(Intercept)"           = "Intercept",
#                   "log_size"              = "log(Size)",
#                   "transitiont2"          = "Transition (t2 vs t1)",
#                   "log_size:transitiont2" = "log(Size) × Transition"),
# 
#     across(c(estimate, std.error, statistic), ~ round(., 3)),
#     p.value = pvalue(p.value, accuracy = 0.001)
#   ) %>%
#   flextable() %>%
#   set_caption("Table X. Logistic regression of coral survival as a function of size and transition") %>%
#   set_header_labels(
#     term       = "Term",
#     estimate   = "Estimate",
#     std.error  = "Std. Error",
#     statistic  = "z value",
#     p.value    = "p-value"
#   ) %>%
#   autofit()
# 
# doc <- read_docx() %>%
#   body_add_par("Survival GLM (s_global)", style = "heading 1") %>%
#   body_add_flextable(s_table)
# 
# # 4. Save Word document
# print(doc, target = "tables/global_survival_table.docx")



### Reproduction ------------------------------------------------------------


r_global  <- glm(repro ~ log_size * transition, family = binomial, data = df_porites)

summary(r_global)

plot(ggpredict(r_global, terms = c("log_size[all]", "transition")))

pred_repro <- ggpredict(r_global, terms = c("log_size[all]", "transition"))
pred_repro_df <- as.data.frame(pred_repro)

pred_repro_df$group <- dplyr::recode(pred_repro_df$group,
                                     "t1" = "2021–2022",
                                     "t2" = "2022–2023")

repro = ggplot(pred_repro_df,
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
  labs(x     = "Size at time t",
       y     = "Probability of reproduction",
       colour = "Transition",
       fill   = "Transition") +
  theme_classic()+
  theme(axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text  = element_text(size = 14),
        legend.text = element_text(size = 12))

# ggsave('figures/results/repro_model.png', units = 'in', width = 10, height = 8, dpi = 600)

# r_tidy <- tidy(r_global)
# 
# # 2. Format and label
# r_table <- r_tidy %>%
#   mutate(
#     term = recode(
#       term,
#       "(Intercept)"           = "Intercept",
#       "log_size"              = "log(Size)",
#       "transitiont2"          = "Transition (t2 vs t1)",
#       "log_size:transitiont2" = "log(Size) × Transition"
#     ),
#     across(c(estimate, std.error, statistic), ~ round(., 3)),
#     p.value = pvalue(p.value, accuracy = 0.001)
#   ) %>%
#   flextable() %>%
#   set_caption("Table X. Logistic regression of reproduction as a function of size and transition") %>%
#   set_header_labels(
#     term      = "Term",
#     estimate  = "Estimate",
#     std.error = "Std. Error",
#     statistic = "z value",
#     p.value   = "p-value"
#   ) %>%
#   autofit()
# 
# doc <- read_docx() %>%
#   body_add_par("Reproduction GLM (r_global)", style = "heading 1") %>%
#   body_add_flextable(r_table)
# 
# print(doc, target = "tables/r_global_repro_table.docx")



### Seeds -------------------------------------------------------------------


seed_global  <- glmmTMB(flower_n_int ~ log_size * transition, family = nbinom2, data = df_porites)
summary(seed_global)
plot(ggpredict(seed_global, terms = c("log_size[all]", "transition")))

pred_seed <- ggpredict(seed_global, terms = c("log_size[all]", "transition"))
pred_seed_df <- as.data.frame(pred_seed)
pred_seed_df$group <- dplyr::recode(pred_seed_df$group,
                                    "t1" = "2021–2022",
                                    "t2" = "2022–2023")

recruit = ggplot(pred_seed_df,
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
  labs(x     = "Size at time t",
       y     = "Expected recruit count",
       colour = "Transition",
       fill   = "Transition") +
  theme_classic()+
  theme(axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text  = element_text(size = 14),
        legend.text = element_text(size = 12))

#ggsave('figures/results/recruit_model.png', units = 'in', width = 10, height = 8, dpi = 600)
# recr_tidy <- tidy(seed_global, effects = "fixed")
# 
# # 2. Format table
# recr_table <- recr_tidy %>%
#   mutate(
#     term = recode(
#       term,
#       "(Intercept)"           = "Intercept",
#       "log_size"              = "log(Size)",
#       "transitiont2"          = "Transition (t2 vs t1)",
#       "log_size:transitiont2" = "log(Size) × Transition"
#     ),
#     across(c(estimate, std.error, statistic), ~ round(., 3)),
#     p.value = pvalue(p.value, accuracy = 0.001)
#   ) %>%
#   flextable() %>%
#   set_caption("Table X. Negative binomial model of recruit counts as a function of colony size and transition") %>%
#   set_header_labels(
#     term      = "Term",
#     estimate  = "Estimate",
#     std.error = "Std. Error",
#     statistic = "z value",
#     p.value   = "p-value"
#   ) %>%
#   autofit()
# 
# # 3. Add to Word document
# doc <- read_docx() %>%
#   body_add_par("Recruitment Model (recruits_global)", style = "heading 1") %>%
#   body_add_flextable(recr_table)
# 
# # 4. Save the Word document
# print(doc, target = "tables/recruits_global_table.docx")



### Plot --------------------------------------------------------------------

ggpubr::ggarrange(growth, 
                  survival,
                  repro, recruit, 
                  ncol = 2, nrow = 4, 
                  labels = c("A)", "B)", "C)", "D)"),
                  common.legend = TRUE, legend = "bottom")

#ggsave('figures/results/global_models_bothv4.png', units = 'in', width = 16, height = 16, dpi = 600)



## IPM ---------------------------------------------------------------------


### Parameters --------------------------------------------------------------

recr_data <- subset(df_porites, action == "born" & area2 < 5)

recr_mu  <- mean(recr_data$log_size_next)
recr_sd  <- sd(recr_data$log_size_next)

grow_sd  <- sd(resid(g_global2))

params_global1 <- list(
  grow_mod   = use_vr_model(g_global2),      
  surv_mod   = s_global,      
  repr_mod   = r_global,    
  seed_mod   = use_vr_model(seed_global),      
  grow_sd    = grow_sd,
  recr_mu    = recr_mu,
  recr_sd    = recr_sd,
  transition = "t1"
)

params_global2 <- list(
  grow_mod   = use_vr_model(g_global2),      # lm with log_size * transition
  surv_mod   = s_global,      # glm with log_size * transition
  repr_mod   = r_global,      # glm binomial with log_size * transition
  seed_mod   = use_vr_model(seed_global),      # glm.nb with log_size * transition
  grow_sd    = grow_sd,
  recr_mu    = recr_mu,
  recr_sd    = recr_sd,
  transition = "t2"
)



### Mesh Size ---------------------------------------------------------------


L <- min(c(df_porites$log_size, df_porites$log_size_next), na.rm = TRUE) * 1.2
U <- max(c(df_porites$log_size, df_porites$log_size_next), na.rm = TRUE) * 1.2

build_ipm_global_mesh <- function(mesh_size, transition, params_global1,
                                  surv_mod, grow_mod, repr_mod, seed_mod,
                                  grow_sd, recr_mu, recr_sd) {
  
  init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "det") %>%
    
    # --- P Kernel (Survival * Growth) ---
    define_kernel(
      name    = "P",
      family  = "CC",
      formula = s * g,
      
      s = predict(surv_mod,
                  newdata = data.frame(
                    log_size   = sa_1,
                    transition = rep(transition, length(sa_1))
                  ),
                  type = "response"),
      
      g_mu = predict(grow_mod,
                     newdata = data.frame(
                       log_size   = sa_1,
                       transition = rep(transition, length(sa_1))
                     ),
                     type = "response"),
      
      g = dnorm(sa_2, g_mu, grow_sd),
      
      states        = list(c("sa")),
      data_list     = params_global1,
      uses_par_sets = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions(fun = "norm",
                                              target = "g")
    ) %>%
    
    # --- F Kernel (Fecundity) ---
    define_kernel(
      name    = "F",
      family  = "CC",
      formula = f_r * f_s * f_d,
      
      f_r = predict(repr_mod,
                    newdata = data.frame(
                      log_size   = sa_1,
                      transition = rep(transition, length(sa_1))
                    ),
                    type = "response"),
      
      f_s = predict(seed_mod,
                    newdata = data.frame(
                      log_size   = sa_1,
                      transition = rep(transition, length(sa_1))
                    ),
                    type = "response"),
      
      f_d = dnorm(sa_2, recr_mu, recr_sd),
      
      states        = list(c("sa")),
      data_list     = params_global1,
      uses_par_sets = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions(fun = "norm",
                                              target = "f_d")
    ) %>%
    
    define_impl(
      make_impl_args_list(
        kernel_names = c("P", "F"),
        int_rule     = rep("midpoint", 2),
        state_start  = rep("sa", 2),
        state_end    = rep("sa", 2)
      )
    ) %>%
    
    define_domains(
      sa = c(L, U, mesh_size)   # <--- mesh_size goes here
    ) %>%
    
    define_pop_state(
      n_sa = runif(mesh_size)   # <--- and here
    ) %>%
    
    make_ipm(iterate = TRUE, iterations = 1000)
}

mesh_seq <- seq(10, 500, by = 20)

lambda_mesh <- map_dfr(mesh_seq, ~{
  ipm_tmp <- build_ipm_global_mesh(
    mesh_size   = .x,
    transition  = transition,   # or "t1"/"t2" if fixed
    params_global1 = params_global1,
    surv_mod    = surv_mod,
    grow_mod    = grow_mod,
    repr_mod    = repr_mod,
    seed_mod    = seed_mod,
    grow_sd     = grow_sd,
    recr_mu     = recr_mu,
    recr_sd     = recr_sd
  )
  
  data.frame(
    mesh_size = .x,
    lambda    = lambda(ipm_tmp)
  )
})

lambda_mesh

ggplot(lambda_mesh, aes(mesh_size, lambda)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(
    name   = "Mesh size",
    limits = c(0, 500),
  ) +
  scale_y_continuous(
    name   = expression(lambda) ) +
  scale_fill_viridis_c(option = "H") +
  theme_classic()+
  theme(axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text  = element_text(size = 14),
        legend.text = element_text(size = 12)) 



### IPM Function ------------------------------------------------------------

L <- min(c(df$log_size, df$log_size_next), na.rm = TRUE) * 1.2
U <- max(c(df$log_size, df$log_size_next), na.rm = TRUE) * 1.2



#### Transition 1 ------------------------------------------------------------


ipm_global1 <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "det") %>%
  
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
    data_list     = params_global1,
    uses_par_sets  = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions(fun = "norm",
                                            target = "g")
  ) %>%
  
  define_kernel(
    name    = "F",
    family  = "CC",
    formula = f_r * f_s * f_d,
    
    f_r = predict(repr_mod,
                  newdata = data.frame(
                    log_size   = sa_1,
                    transition = rep(transition, length(sa_1))
                  ),
                  type = "response"),
    
    f_s = predict(seed_mod,
                  newdata = data.frame(
                    log_size   = sa_1,
                    transition = rep(transition, length(sa_1))
                  ),
                  type = "response"),
    
    f_d = dnorm(sa_2, recr_mu, recr_sd),
    
    
    states        = list(c("sa")),
    data_list     = params_global1,
    uses_par_sets  = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions(fun = "norm",
                                            target = "f_d")
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P", "F"),
      int_rule     = rep('midpoint', 2),
      state_start    = rep("sa", 2),
      state_end      = rep("sa", 2)
    )
  ) %>%
  define_domains(
    sa = c(L,
           U,
           200)
  ) %>%
  define_pop_state(n_sa = runif(200)) %>%
  make_ipm()


#### Transition 2 ------------------------------------------------------------


ipm_global2 <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "det") %>%
  
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
    data_list     = params_global2,
    uses_par_sets  = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions(fun = "norm",
                                            target = "g")
  ) %>%
  
  define_kernel(
    name    = "F",
    family  = "CC",
    formula = f_r * f_s * f_d,
    
    f_r = predict(repr_mod,
                  newdata = data.frame(
                    log_size   = sa_1,
                    transition = rep(transition, length(sa_1))
                  ),
                  type = "response"),
    
    f_s = predict(seed_mod,
                  newdata = data.frame(
                    log_size   = sa_1,
                    transition = rep(transition, length(sa_1))
                  ),
                  type = "response"),
    
    f_d = dnorm(sa_2, recr_mu, recr_sd),
    
    
    states        = list(c("sa")),
    data_list     = params_global2,
    uses_par_sets  = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions(fun = "norm",
                                            target = "f_d")
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P", "F"),
      int_rule     = rep('midpoint', 2),
      state_start    = rep("sa", 2),
      state_end      = rep("sa", 2)
    )
  ) %>%
  define_domains(
    sa = c(L2,
           U2,
           200)
  ) %>%
  define_pop_state(n_sa = runif(200)) %>%
  make_ipm()


#### Lambda ------------------------------------------------------------------



lambda(ipm_global1)
lambda(ipm_global2)


## Lambda Bootstrap ---------------------------------------------------------------

# ── Bootstrap function for Porites (P + F kernels) ────────────────────────────
# build_ipm_lambda_por <- function(df_boot, transition = "t1", mesh_size = 200) {
#   
#   # Refit vital rate models
#   g_mod_boot    <- update(g_global2,   data = df_boot)
#   s_mod_boot    <- update(s_global,    data = df_boot)
#   r_mod_boot    <- update(r_global,    data = df_boot)
#   seed_mod_boot <- update(seed_global, data = df_boot)
#   
#   grow_sd_boot <- sd(resid(g_mod_boot), na.rm = TRUE)
#   
#   # Recruitment parameters from bootstrap sample
#   recr_data_boot <- subset(df_boot, action == "born" & area2 < 5)
#   recr_mu_boot   <- mean(recr_data_boot$log_size_next, na.rm = TRUE)
#   recr_sd_boot   <- sd(recr_data_boot$log_size_next,   na.rm = TRUE)
#   
#   params_boot <- list(
#     grow_mod   = use_vr_model(g_mod_boot),
#     surv_mod   = s_mod_boot,
#     repr_mod   = r_mod_boot,
#     seed_mod   = use_vr_model(seed_mod_boot),
#     grow_sd    = grow_sd_boot,
#     recr_mu    = recr_mu_boot,
#     recr_sd    = recr_sd_boot,
#     transition = transition
#   )
#   
#   # Shared domain from full dataset
#   L_boot <- min(c(df_porites$log_size, df_porites$log_size_next), na.rm = TRUE) * 1.2
#   U_boot <- max(c(df_porites$log_size, df_porites$log_size_next), na.rm = TRUE) * 1.2
#   
#   ipm_boot <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "det") %>%
#     define_kernel(
#       name = "P", family = "CC", formula = s * g,
#       s    = predict(surv_mod,
#                      newdata = data.frame(log_size   = sa_1,
#                                           transition = rep(transition, length(sa_1))),
#                      type = "response"),
#       g_mu = predict(grow_mod,
#                      newdata = data.frame(log_size   = sa_1,
#                                           transition = rep(transition, length(sa_1))),
#                      type = "response"),
#       g    = dnorm(sa_2, g_mu, grow_sd),
#       states        = list(c("sa")),
#       data_list     = params_boot,
#       uses_par_sets = FALSE,
#       evict_cor     = TRUE,
#       evict_fun     = truncated_distributions(fun = "norm", target = "g")
#     ) %>%
#     define_kernel(
#       name = "F", family = "CC", formula = f_r * f_s * f_d,
#       f_r  = predict(repr_mod,
#                      newdata = data.frame(log_size   = sa_1,
#                                           transition = rep(transition, length(sa_1))),
#                      type = "response"),
#       f_s  = predict(seed_mod,
#                      newdata = data.frame(log_size   = sa_1,
#                                           transition = rep(transition, length(sa_1))),
#                      type = "response"),
#       f_d  = dnorm(sa_2, recr_mu, recr_sd),
#       states        = list(c("sa")),
#       data_list     = params_boot,
#       uses_par_sets = FALSE,
#       evict_cor     = TRUE,
#       evict_fun     = truncated_distributions(fun = "norm", target = "f_d")
#     ) %>%
#     define_impl(make_impl_args_list(
#       kernel_names = c("P", "F"),
#       int_rule     = rep("midpoint", 2),
#       state_start  = rep("sa", 2),
#       state_end    = rep("sa", 2)
#     )) %>%
#     define_domains(sa = c(L_boot, U_boot, mesh_size)) %>%
#     define_pop_state(n_sa = runif(mesh_size)) %>%
#     make_ipm(iterate = TRUE, iterations = 1000)
#   
#   lambda(ipm_boot)
# }
# 
# # ── One bootstrap replicate returning both t1 and t2 ──────────────────────────
# one_boot_por <- function(i) {
#   set.seed(i)
#   idx     <- sample(seq_len(nrow(df_porites)), replace = TRUE)
#   df_boot <- df_porites[idx, ]
#   
#   # Guard against bootstrap samples with no recruits
#   if (nrow(subset(df_boot, action == "born" & area2 < 5)) < 2) {
#     return(c(t1 = NA_real_, t2 = NA_real_))
#   }
#   
#   lam_t1 <- tryCatch(suppressWarnings(build_ipm_lambda_por(df_boot, transition = "t1")),
#                      error = function(e) NA_real_)
#   lam_t2 <- tryCatch(suppressWarnings(build_ipm_lambda_por(df_boot, transition = "t2")),
#                      error = function(e) NA_real_)
#   
#   c(t1 = lam_t1, t2 = lam_t2)
# }
# 
# # ── Windows parallel setup ─────────────────────────────────────────────────────
# n_cores <- detectCores() - 6
# cl      <- makeCluster(n_cores)
# 
# clusterExport(cl, varlist = c("df_porites",
#                               "build_ipm_lambda_por",
#                               "g_global2",
#                               "s_global",
#                               "r_global",
#                               "seed_global",
#                               "use_vr_model"))
# 
# clusterEvalQ(cl, {
#   library(ipmr)
#   library(dplyr)
#   library(glmmTMB)
#   
# })
# 
# # ── Run with progress bar ──────────────────────────────────────────────────────
# set.seed(123)
# cat("Starting bootstrap — running 1200 replicates across", n_cores, "cores\n")
# t_start <- proc.time()
# 
# results <- pblapply(1:1200, one_boot_por, cl = cl)
# 
# t_end <- proc.time()
# stopCluster(cl)
# cat("Done! Elapsed time:", round((t_end - t_start)["elapsed"] / 60, 1), "minutes\n")
# 
# # ── Parse results ──────────────────────────────────────────────────────────────
# results_df <- do.call(rbind, lapply(results, function(x) {
#   if (is.null(x) || length(x) != 2) return(c(t1 = NA_real_, t2 = NA_real_))
#   x
# }))
# 
# colnames(results_df) <- c("t1", "t2")
# 
# cat("Valid t1 replicates:", sum(is.finite(results_df[, "t1"])), "\n")
# cat("Valid t2 replicates:", sum(is.finite(results_df[, "t2"])), "\n")
# 
# lambda_boot_t1 <- na.omit(results_df[, "t1"])[1:1000]
# lambda_boot_t2 <- na.omit(results_df[, "t2"])[1:1000]
# 
# # ── Summaries ─────────────────────────────────────────────────────────────────
# ci_t1 <- quantile(lambda_boot_t1, c(0.025, 0.5, 0.975), na.rm = TRUE)
# ci_t2 <- quantile(lambda_boot_t2, c(0.025, 0.5, 0.975), na.rm = TRUE)
# 
# mean_t1 <- mean(lambda_boot_t1, na.rm = TRUE)
# mean_t2 <- mean(lambda_boot_t2, na.rm = TRUE)
# mean_t1; ci_t1
# mean_t2; ci_t2
# 
# # ── Histograms ────────────────────────────────────────────────────────────────
# par(mfrow = c(1, 2))
# hist(lambda_boot_t1, breaks = 30,
#      main = expression("Bootstrap " * lambda * " 2021-2022"),
#      xlab = expression(lambda))
# abline(v = ci_t1[c(1,3)], col = "red",  lty = 2)
# abline(v = ci_t1[2],       col = "blue", lwd = 2)
# 
# hist(lambda_boot_t2, breaks = 30,
#      main = expression("Bootstrap " * lambda * " 2022-2023"),
#      xlab = expression(lambda),
#      xlim = c(min(lambda_boot_t2), max(lambda_boot_t2) * 1.01))
# abline(v = ci_t2[c(1,3)], col = "red",  lty = 2)
# abline(v = ci_t2[2],       col = "blue", lwd = 2)
# par(mfrow = c(1, 1))
# 
# # ── Combined data frame for ggplot ────────────────────────────────────────────
# all_lambdas_por <- data.frame(
#   lambda     = c(lambda_boot_t1, lambda_boot_t2),
#   transition = rep(c("2021–2022", "2022–2023"), each = 1000)
# )
# 
# 
# write_csv(all_lambdas_por,"Markdowns/data/global_ipm_lambda_bootstrap.csv")

### Boostrap data -----------------------------------------------------------


all_lambdas = read_csv("Markdowns/data/global_ipm_lambda_bootstrap.csv")  

### Density Plot ###

boot = ggplot(all_lambdas, aes(x = lambda, fill = transition)) +
  geom_density(alpha = 0.35) +
  labs(
    x = expression(lambda),
    y = "Density", fill = "Transition") +
  scale_fill_viridis_d(option = "H", begin = .9, end = 0.2) +
  theme_classic()+
  theme(axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text  = element_text(size = 14),
        legend.text = element_text(size = 12))

boot

### Stochastic lambda projection --------------------------------------------


# settings
n_sims  <- 500
start_year <- 2025
n_years <- 51
years <- start_year:(start_year + n_years - 1)
N0 <- 150

boot <- all_lambdas %>%
  dplyr::select(transition, lambda) %>%
  filter(!is.na(lambda))

sim <- boot %>%
  group_by(transition) %>%
  summarise(lambda_draws = list(lambda), .groups = "drop") %>%
  tidyr::uncount(n_sims, .id = "sim_id") %>%
  rowwise() %>%
  mutate(
    lambda_path = list(sample(lambda_draws, size = n_years - 1, replace = FALSE)),
    year = list(years),
    N = list(c(
      N0,
      N0 * exp(cumsum(log(lambda_path)))
    ))
  ) %>%
  unnest(c(year, N)) %>%
  ungroup()

# Plot
f2a = ggplot(sim, aes(x = year, y = N,
                      group = interaction(transition, sim_id),
                      color = transition)) +
  geom_line(alpha = .1, linewidth = 0.3) +
  # annotation_custom(
  #   grob_imgpast,
  #   xmin = 2067, xmax = 2072,   # <-- adjust to your x scale
  #   ymin = 2000,
  #   ymax = 2300
  # )+
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
f2a



## Generation time ---------------------------------------------------------

R_nought<-function(ipm_obj) {
  Pm<-ipm_obj$sub_kernels$P 
  Fm<-ipm_obj$sub_kernels$F 
  I <-diag(dim(Pm)[1]) 
  N <-solve(I-Pm) 
  R <-Fm %*%N 
  return( Re(eigen(R)$values)[1] 
  ) 
}

# Generation Time function: How many years to replace themselves demographically

gen_time<-function(ipm_obj) { 
  lamb <-unname(lambda(ipm_obj)) 
  r_nought <-R_nought(ipm_obj) 
  return(log(r_nought)/log(lamb)) 
}


R_nought(ipm_global1)
R_nought(ipm_global2)
gen_time(ipm_global1)
gen_time(ipm_global2)



## Heat Maps ---------------------------------------------------------------


### P kernel ----------------------------------------------------------------


#### Transition 1 ------------------------------------------------------------


Kt1 <- ipm_global1$sub_kernels$P

wt1 <- Re(eigen(Kt1)$vectors[, 1])          # right eigenvector
vt1 <- Re(eigen(t(Kt1))$vectors[, 1])       # left eigenvector

# Ensure positive orientation
wt1 <- abs(wt1)
vt1 <- abs(vt1)

mesh_info_port1<-int_mesh(ipm_global1) 

d_sat1 <- mesh_info_port1$d_sa

sens_mat_port1 <- outer(vt1, wt1) / sum(vt1 * wt1 * d_sat1)
elas_mat_port1 <- sens_mat_port1 * Kt1 / lambda(ipm_global1)

sens_df_port1 <- ipm_to_df(sens_mat_port1) 
elas_df_port1 <- ipm_to_df(elas_mat_port1)

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


p_df_port1 <- data.frame(
  t     = mesh_info_port1$sa_1,
  t_1   = mesh_info_port1$sa_2,
  value = as.vector(ipm_global1$sub_kernels$P)
)

p_plt_port1 <- ggplot(p_df_port1) + 
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

p_plt_port1

# 1) Conditional *mean* size at t+1 for each t
ridge_mean_port1 <- p_df_port1 |>
  dplyr::group_by(t) |>
  dplyr::summarise(t1_bar = sum(t_1 * value) / sum(value), .groups = "drop") |>
  dplyr::mutate(increment = t1_bar - t)


# 2) Overlay that mean curve
p_plt_port1 +
  geom_line(data = ridge_mean_port1, aes(x = t, y = t1_bar),
            linewidth = 1.1, color = "white")




#### Transition 2 ------------------------------------------------------------


Kt2 <- ipm_global2$sub_kernels$P

wt2 <- Re(eigen(Kt2)$vectors[, 1])          # right eigenvector
vt2 <- Re(eigen(t(Kt2))$vectors[, 1])       # left eigenvector

# Ensure positive orientation
wt2 <- abs(wt2)
vt2 <- abs(vt2)

mesh_info_port2<-int_mesh(ipm_global2) 

d_sat2 <- mesh_info_port2$d_sa

sens_mat_port2 <- outer(vt2, wt2) / sum(vt2 * wt2 * d_sat2)
elas_mat_port2 <- sens_mat_port2 * Kt2 / lambda(ipm_global2)

sens_df_port2 <- ipm_to_df(sens_mat_port2) 
elas_df_port2 <- ipm_to_df(elas_mat_port2)

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


p_df_port2 <- data.frame(
  t     = mesh_info_port2$sa_1,
  t_1   = mesh_info_port2$sa_2,
  value = as.vector(ipm_global2$sub_kernels$P)
)

p_plt_port2 <- ggplot(p_df_port2) + 
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

p_plt_port2

# 1) Conditional *mean* size at t+1 for each t
ridge_mean_port2 <- p_df_port2 |>
  dplyr::group_by(t) |>
  dplyr::summarise(t1_bar = sum(t_1 * value) / sum(value), .groups = "drop") |>
  dplyr::mutate(increment = t1_bar - t)


# 2) Overlay that mean curve
p_plt_port2 +
  geom_line(data = ridge_mean_port2, aes(x = t, y = t1_bar),
            linewidth = 1.1, color = "white")

