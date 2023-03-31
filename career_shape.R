library(tidyverse)
library(psych)      # for PCA
library(factoextra) # for clustering
library(lubridate)  # for dates
library(tidymodels)
library(brms)       # for model building and comparisons
library(ggdist)
library(tidybayes)
library(modelr)
library(cowplot)    # for plot tiling
library(ggrepel)    # for labeling points
library(ggbeeswarm) # for geom_quasirandom
library(afex)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load + Wrangle Data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

setwd('data')
filenames <- list.files()
d <- lapply(filenames, read_csv)

# I accidentally called all the embedding variables almog_embedding_1 etc. Let's just fix that.
remove_almog <- function(x) rename_with(x, ~str_remove(.x, 'almog_'), starts_with("almog_"))
d <- lapply(d, remove_almog)

scholars_ordered <- c("Almog Simchon", "Daniel L Schacter", "Duncan Astle", "Falk Eippert",
                      "Andrea Berger", "Joshua D Greene", "John Duncan", "Lars Meyer",
                      "Florina Uzefovsky", "Samuel J Gershman", "Kate Baker", "Roland G Benoit",
                      "Nachshon Meiran", "Steven Pinker", "Michael Anderson", "Stephanie Theves",
                      "Yoav Kessler", "Tomer D Ullman", "Tim Dalgleish", "Veronika Engert")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Principal Component Analysis and t-SNE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
add_pca = function(dat) {
  pca <- dat %>% select(embedding_1:embedding_768) %>% as.matrix() %>% principal(nfactors = 2)
  scholar_name = paste0(str_replace(str_to_lower(dat$scholar_name[1]), ' ', '_'), "_pca")
  assign(scholar_name, pca, envir = .GlobalEnv)
  dat$pc1 <- pca$scores[,"RC1"]
  dat$pc2 <- pca$scores[,"RC2"]
  dat
}

d <- lapply(d, add_pca)

# Print Proportion Variance Explained
lapply(scholars_ordered[19:20], function(x){print(x)
  get(paste0(str_replace(str_to_lower(x), ' ', '_'), "_pca"))$Vaccounted})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# One Tibble to Rule them All
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
d_agg <- bind_rows(d)

d_agg <- d_agg %>% 
  rename(pub_index = ...1) %>% 
  # Remove outliers from the beginning of the career 
  # (A few of these are clearly mistakes by the Semantic Scholar Algorithm)
  left_join(d_agg %>% 
              group_by(scholar_name) %>% 
              summarise(min_date = quantile(date, probs = c(.01), na.rm = T)))  %>% 
  filter(date > min_date)

# New Dataset with Relevant Vars Aggregated by scholar_name
d_analysis <- d_agg %>% 
  # Remove outliers from the beginning of the career 
  left_join(d_agg %>% 
              group_by(scholar_name) %>% 
              summarise(min_date = quantile(date, probs = c(.01), na.rm = T))) %>% 
  filter(date > min_date) %>% 
  # Aggregated Vars
  group_by(group, scholar_name) %>% 
  summarise(min_date = min(date, na.rm = T),
            max_date = max(date, na.rm = T),
            career_length_yrs = as.numeric(max_date-min_date, units = "days")/365,
            npubs = n(),
            pubfreq = npubs/career_length_yrs,
            # Average Number of Co-Authors
            n_authors = mean(n_authors, na.rm = T)) %>% 
  ungroup() %>% 
  # Gender
  mutate(gender = if_else(scholar_name %in% c("Almog Simchon", "Daniel L Schacter", "Duncan Astle", "Falk Eippert",
                                              "Joshua D Greene", "John Duncan", "Lars Meyer", "Samuel J Gershman", 
                                              "Roland G Benoit", "Nachshon Meiran", "Steven Pinker", "Michael Anderson", 
                                              "Yoav Kessler", "Tomer D Ullman", "Tim Dalgleish"),
                          "M", "F"),
         group_type = if_else(group %in% c("Ben Gurion University Psychology Department", "Harvard CBB Group"),
                              "University", "Institute")) %>% 
  # Optimal K
  left_join(tibble(scholar_name = scholars_ordered,
                   optimal_k = unlist(lapply(gap_stats, function(i){cluster::maxSE(i$Tab[,"gap"], i$Tab[,"SE.sim"])}))
  )
  ) %>% 
  # Euclidean Distance Between Average of 1st and last 3rd of Career
  left_join(
    d_agg %>% 
      group_by(scholar_name) %>% 
      mutate(third = as.numeric(cut_interval(date, 3))) %>% 
      group_by(scholar_name, third) %>% 
      summarise_at(vars(embedding_1:embedding_768), mean) %>%
      filter(third != 2) %>% 
      select(scholar_name, embedding_1:embedding_768) %>% 
      group_by(scholar_name) %>% 
      summarise_at(vars(embedding_1:embedding_768), diff) %>% 
      ungroup() %>%
      mutate_at(vars(embedding_1:embedding_768), ~.^2) %>% 
      mutate(euc_dist = sqrt(rowSums(select(., -scholar_name))))
  ) %>% 
  # Total Variance (i.e. sum of Eigenvalues in PCA)
  left_join(
    d_agg %>% 
      group_by(scholar_name) %>% 
      summarise_at(vars(embedding_1:embedding_768), var, na.rm = TRUE) %>% 
      mutate(total_var = rowSums(select(., embedding_1:embedding_768))) %>% 
      select(scholar_name, total_var)
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Exploratory Graphical Analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

## Pretty Plot of PC1 by PC2 for each Scholar, Colored by Years Since First Publication
p1 <- d_agg %>% 
  # Coloring Variable: Years Since first Publication
  group_by(scholar_name) %>% 
  mutate(min_date = min(date, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(career_length = as.numeric(date - min_date, units="days")) %>% 
  # Count Publications - for Subtly Varying Opacity as Needed
  left_join(d_agg %>% group_by(scholar_name) %>% summarise(n = n())) %>% 
  #Plot
  ggplot(aes(pc1, pc2, color = career_length/365, alpha = 1/log(n))) +
    geom_point() +
    facet_wrap(~ fct_reorder(scholar_name, as.numeric(as_factor(group))),
               ncol = 4, dir = "v") +
    scale_color_gradient2(low = "blue", mid = "red", high = "orange", midpoint = 20) + 
    theme_bw() +
    guides(alpha = "none") +
    labs(title = "", x = "Subject Matter (PC1)", y = "Subject Matter (PC2)", color = "Years Since\nFirst Publication") + 
    theme(axis.text = element_blank(), strip.background = element_rect(fill = "lavender")) 
    
  # Annotate Groups (Hacky Solution)
  group_labels <- tibble(group = c("Ben Gurion University", "Harvard University", 
                   "CBSU", "Max Planck Institute"),
                   color = c(0, 10, 30, 40)) %>% 
    ggplot() + 
    geom_text(aes(x = 0, y = 0, label = color, color = color)) + 
    geom_text(aes(x = 0, y = 0, label = group), color = "black", fontface = "bold", size = 4.5) + 
    scale_color_gradient(low = "white", high = "white") +
    facet_wrap(~ group, nrow = 1) + theme_void() +
    theme(strip.text = element_blank(),
          legend.text = element_text(color = "white"),
          legend.title = element_text(color = "white"),
          legend.key = element_rect(fill = "white", linewidth = 0))
  
  # Final Plot       
  plot_grid(group_labels, p1, nrow = 2, rel_heights = c(1, 25), align = 'v')

# Plot Careers Alongside Co-Authors
plot_grid(group_labels, 
  d_agg %>% 
  ggplot(aes(date, group = scholar_name)) +
      geom_smooth(aes(y = n_authors-1, color = "Number of Co-Authors"), 
                  method = "loess", method.args = list(degree = 1), se = F) + 
      # geom_line(aes(y = n_citations, color = "n_citations"), se = F) + 
      geom_smooth(aes(y = pc1/8, color = "PC1 (rescaled)"), 
                  method = "loess", method.args = list(degree = 1), se = F) + 
      facet_wrap(~ fct_reorder(scholar_name, as.numeric(as_factor(group))),
                 ncol = 4, dir = "v") +
      labs(x = "", y = "", color = "", caption = "Lines are smoothed using locally estimated scatterplot smoothing (LOESS) with degree 1 polynomials.") +
      coord_cartesian(ylim=c(-15, 40)) +
      theme_bw() +
      theme(strip.background = element_rect(fill = "lavender")),
  nrow = 2, rel_heights = c(1, 25), align = 'v')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# K-Means Clustering
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

find_optimal_k = function(name) {
  cluster::clusGap(d_agg %>% filter(scholar_name == name) %>% select(embedding_1:embedding_768),
                   kmeans, K.max = 10, B = 500)
}

gap_stats <- lapply(scholars_ordered, find_optimal_k)

# Visualize Optimal Clustering
plot_optimal_k = function(gap_stat) {
  fviz_gap_stat(gap_stat) + labs(x = "k", y = "", title = "") + theme(axis.text.y = element_blank())
}

plotlist <- lapply(scholars_ordered, plot_optimal_k)

plot_grid(tibble(group = c("Ben Gurion University", "Harvard University", 
                           "CBSU", "Max Planck Institute")) %>% 
            ggplot() + geom_text(aes(x = 0, y = 0, label = group), color = "black", fontface = "bold", size = 4.5) + 
            facet_wrap(~ group, nrow = 1) + theme_void() + theme(strip.text = element_blank()), 
          plot_grid(plotlist = plotlist, labels = scholars_ordered, label_size = 12, label_fontface = "plain", ncol = 4), 
          ncol = 1, rel_heights = c(2, 25))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Reorganized Dataset for Modeling Over Time
d_threeyear <- d_agg %>% 
  group_by(group, scholar_name) %>% 
  # Three year segments
  mutate(threeyear = as.integer(cut_width(date(date), 1095))) %>% 
  group_by(group, scholar_name, threeyear) %>% 
  summarise(across(embedding_1:embedding_768, ~mean(.x, na.rm = TRUE), .names = "{.col}_diff"),
            across(embedding_1:embedding_768, ~var(.x, na.rm = TRUE), .names = "{.col}_var"),
            n_pubs = n(),
            n_authors = sum((n_authors-1), na.rm = TRUE)/n()) %>% 
  mutate(across(embedding_1_diff:embedding_768_diff, ~diff(c(NA, .x)))) %>% 
  mutate(across(embedding_1_diff:embedding_768_diff, ~.x^2)) %>% 
  rowwise() %>% 
  mutate(euc_dist = sqrt(sum(across(embedding_1_diff:embedding_768_diff), na.rm = TRUE)),
         total_var = sum(across(embedding_1_var:embedding_768_var), na.rm = TRUE)) %>% 
  # Add Gender
  left_join(d_analysis %>% select(scholar_name, gender, group_type)) %>% 
  select(group, scholar_name, gender, group_type, total_var, euc_dist, n_pubs, n_authors, threeyear) %>% 
  # remove 1st (always 0 difference) and last row (likely smaller n)
  group_by(scholar_name) %>% 
  slice(2:(n()-1))

# Bayesian Multilevel Polynomial Model Comparisons
reduced_mod <- brm(euc_dist ~ threeyear + (threeyear | scholar_name), data = d_threeyear, 
                   prior = c(
                     prior(normal(0, 50), class = "b", coef = "threeyear")
                   ),
                   chains = 4, cores = 4)
quadratic_mod <- brm(euc_dist ~ poly(threeyear,2) + (poly(threeyear,2) | scholar_name), data = d_threeyear, 
                     prior = c(
                       prior(normal(0, 50), class = "b", coef = "polythreeyear21"),
                       prior(normal(0, 50), class = "b", coef = "polythreeyear22")
                     ),
                     chains = 4, cores = 4)
cubic_mod <- brm(euc_dist ~ poly(threeyear,3) + (poly(threeyear,3) | scholar_name), data = d_threeyear, 
                 prior = c(
                   prior(normal(0, 50), class = "b", coef = "polythreeyear31"),
                   prior(normal(0, 50), class = "b", coef = "polythreeyear32"),
                   prior(normal(0, 50), class = "b", coef = "polythreeyear33")
                 ),
                 chains = 4, cores = 4, iter = 10000, control = list(adapt_delta = .98))
quartic_mod <- brm(euc_dist ~ poly(threeyear,4) + (poly(threeyear,4) | scholar_name), data = d_threeyear, 
                   prior = c(
                     prior(normal(0, 50), class = "b", coef = "polythreeyear41"),
                     prior(normal(0, 50), class = "b", coef = "polythreeyear42"),
                     prior(normal(0, 50), class = "b", coef = "polythreeyear43"),
                     prior(normal(0, 50), class = "b", coef = "polythreeyear44")
                   ),
                   chains = 4, cores = 4, iter = 10000, control = list(adapt_delta = .99))
quintic_mod <- brm(euc_dist ~ poly(threeyear,5) + (poly(threeyear,5) | scholar_name), data = d_threeyear, 
                   prior = c(
                     prior(normal(0, 50), class = "b", coef = "polythreeyear51"),
                     prior(normal(0, 50), class = "b", coef = "polythreeyear52"),
                     prior(normal(0, 50), class = "b", coef = "polythreeyear53"),
                     prior(normal(0, 50), class = "b", coef = "polythreeyear54"),
                     prior(normal(0, 50), class = "b", coef = "polythreeyear55")
                   ),
                   chains = 4, cores = 4, iter = 10000, control = list(adapt_delta = .99))

# Add 10-fold cross-validation
reduced_mod <- add_criterion(reduced_mod, "kfold")
quadratic_mod <- add_criterion(quadratic_mod, "kfold")
cubic_mod <- add_criterion(cubic_mod, "kfold")
quartic_mod <- add_criterion(quartic_mod, "kfold")
quintic_mod <- add_criterion(quintic_mod, "kfold")

pairs(quintic_mod_gender)

loo_compare(reduced_mod, quadratic_mod, cubic_mod, quartic_mod, quintic_mod, criterion = "kfold") # Quartic Wins!

# Try adding n_authors
quartic_mod_coauthors <- brm(euc_dist ~ n_authors + poly(threeyear,4) + (poly(threeyear,4) | scholar_name), data = d_threeyear, 
                   prior = c(
                     prior(normal(0, 50), class = "b")
                   ),
                   chains = 4, cores = 4, iter = 10000, control = list(adapt_delta = .99))

quartic_mod_coauthors <- add_criterion(quartic_mod_coauthors, "kfold")
loo_compare(quartic_mod, quartic_mod_coauthors, criterion = "kfold") # does not improve the model’s expected predictive accuracy 


# Plot euc_dist with respect to time, with posterior predictive of quartic_mod 
d_threeyear %>%
  mutate(gender = if_else(gender == "M", "Men", "Women")) %>% 
  group_by(scholar_name, gender, group_type) %>%
  data_grid(threeyear = seq_range(threeyear, n = 50)) %>%
  add_epred_draws(quartic_mod) %>%
  ggplot(aes(threeyear, euc_dist, color = fct_reorder(scholar_name, as.numeric(as_factor(gender))), group = scholar_name, fill = fct_reorder(scholar_name, as.numeric(as_factor(gender))))) +
  stat_lineribbon(aes(y = .epred), alpha = .6, .width = c(.95)) +
  geom_line(alpha = .8, data = d_threeyear %>% mutate(gender = if_else(gender == "M", "Men", "Women"))) +
  geom_point(data = d_threeyear %>% mutate(gender = if_else(gender == "M", "Men", "Women"))) +
  labs(x = "Years Since First Publication", 
       y = str_wrap("Change in Average Subject Matter Relative to Previous Three Years (Euclidean Distance)", 45),
       caption = str_wrap("Smooth Lines and 95% credible intervals represent expectations of the posterior predictive of a fourth-degree polynomial regression fit, with all parameters allowed to vary by individual.", 93)) +
  theme_bw() +
  facet_grid(vars(group_type), vars(gender)) +
  scale_fill_manual(guide = "none",
                    values = colorRampPalette(RColorBrewer::brewer.pal(8, "Blues"))(20)) +
  # scale_fill_brewer(name = "Confidence", palette = "Greys") +
  scale_color_manual(guide = "none",
                     values = colorRampPalette(RColorBrewer::brewer.pal(8, "Blues"))(20)) +
  scale_x_continuous(breaks = 2:15,
                     labels = seq(6, 45, 3)) +
  theme(strip.background = element_rect(fill = "lavender"))

# Similar process, this time predicting publication frequency by time (using frequentist methods for now)
reduced_mod_freq_npubs <- lme4::lmer(n_pubs ~ threeyear + (threeyear | scholar_name), data = d_threeyear)
quadratic_mod_freq_npubs <- lme4::lmer(n_pubs ~ poly(threeyear,2) + (poly(threeyear,2) | scholar_name), data = d_threeyear)
cubic_mod_freq_npubs <- lme4::lmer(n_pubs ~ poly(threeyear,3) + (poly(threeyear,3) | scholar_name), data = d_threeyear)
quartic_mod_freq_npubs <- lme4::lmer(n_pubs ~ poly(threeyear,4) + (poly(threeyear,4) | scholar_name), data = d_threeyear)

BIC(reduced_mod_freq_npubs, 
    quadratic_mod_freq_npubs,
    cubic_mod_freq_npubs,
    quartic_mod_freq_npubs) # Quadratic Wins!

# Sister plot for above, this time with publication frequency on the y axis
d_threeyear %>%
  mutate(gender = if_else(gender == "M", "Men", "Women")) %>% 
  group_by(scholar_name, gender, group_type) %>%
  data_grid(threeyear = seq_range(threeyear, n = 50)) %>%
  add_predictions(quadratic_mod_freq_npubs) %>% 
  ggplot(aes(threeyear, n_pubs, color = scholar_name, group = scholar_name, fill = fct_reorder(scholar_name, as.numeric(as_factor(gender))))) +
  geom_line(alpha = .8, data = d_threeyear %>% mutate(gender = if_else(gender == "M", "Men", "Women"))) +
  geom_point(data = d_threeyear %>% mutate(gender = if_else(gender == "M", "Men", "Women"))) +
  geom_line(aes(y = pred), linewidth = 1) +
  labs(x = "Years Since First Publication", 
       y = str_wrap("Number of Publications (Three Year Sum)", 45),
       caption = str_wrap("Smooth Lines represent predictions of a second-degree polynomial regression fit, with all parameters allowed to vary by individual.", 70)) +
  theme_bw() +
  scale_fill_manual(guide = "none",
                    values = rev(colorRampPalette(RColorBrewer::brewer.pal(8, "Blues"))(20))) +
  # scale_fill_brewer(name = "Confidence", palette = "Greys") +
  scale_color_manual(guide = "none",
                     values = rev(colorRampPalette(RColorBrewer::brewer.pal(8, "Blues"))(20))) +
  scale_x_continuous(breaks = 2:15,
                     labels = seq(6, 45, 3)) +
  theme(strip.background = element_rect(fill = "lavender"))

d_threeyear %>% filter(n_pubs < 2)
  # The truth comes out! The three authors with the strange shape are also the only three with bins containing only one publication!


# Who is More Likely to Change Interests over their Career?
  # Men or Women, University or Institute
  d_analysis %>% 
    group_by(group_type, gender) %>% 
    summarise(euc_dist = mean(euc_dist)) %>% 
    ggplot(aes(group_type, euc_dist, group = gender, color = gender)) +
    geom_line(size = 1.5) +
    geom_quasirandom(width = .1, data = d_analysis, dodge.width=.4, alpha = .8) + 
    geom_text_repel(aes(label = scholar_name), 
                    size = 3,
                    data = d_analysis) + 
    theme_minimal() +
    theme(axis.text.y = element_blank()) +
    labs(x = "", y = "Semantic Distance Between First and Last Third of Career", color = "Gender")
  
  intrest_change_anova <- aov_ez(id = "scholar_name", dv = "euc_dist",
                                 d_analysis, between = c("gender", "group_type"))
  emmeans::emmeans(intrest_change_anova, c("gender", "group_type"))
  
  # Average Number of Co-Authors (check for interaction with above)
    # In the future, find a way to distinguish student co-authors from colleagues in the field
  d_analysis %>% 
    ggplot(aes(n_authors-1, euc_dist, color = group_type)) +
    geom_point(aes(size = career_length_yrs), alpha = .7) +
    geom_text_repel(aes(label = if_else(euc_dist > 27 | n_authors > 6, scholar_name, "")), 
                    color = "black", size = 3.2, fontface = "italic",
                    data = d_analysis) +
    theme_minimal() +
    scale_color_brewer(palette = "Dark2") +
    theme(axis.text.y = element_blank()) +
    labs(x = "Average Number of Co-Authors", 
         y = "Semantic Distance Between First and Last Third of Career", 
         color = "", size = "Career Length\n(Years)")

  model.comparison(lm(euc_dist ~ n_authors + group_type, d_analysis),
                   lm(euc_dist ~ n_authors*group_type, d_analysis))
  
  cor(d_analysis$n_authors, d_analysis$euc_dist)
# Who is More Likely to Have Multiple Research Interests?
  # By optimal K
    # Men or Women, University or Institute
    d_analysis %>% 
      group_by(group_type, gender) %>% 
      summarise(optimal_k = mean(optimal_k)) %>% 
      ggplot(aes(group_type, optimal_k, group = gender, color = gender)) +
      geom_line(size = 1.5) +
      geom_quasirandom(width = .1, data = d_analysis, dodge.width=.4, alpha = .8) + 
      geom_text_repel(aes(label = scholar_name),
                      size = 3,
                      data = d_analysis) +
      theme_minimal() +
      theme(axis.text.y = element_blank()) +
      labs(x = "", y = "Number of Research Topics (Optimal K)", color = "Gender")
    
    optimal_k_anova <- aov_ez(id = "scholar_name", dv = "optimal_k",
                                   d_analysis, between = c("gender", "group_type"),
                              anova_table = list(es = "pes"))
    optimal_k_anova
    emmeans::emmeans(optimal_k_anova, c("gender", "group_type"))
  
    # Average Number of Co-Authors (check for interaction with above)
    d_analysis %>% 
      ggplot(aes(n_authors-1, optimal_k, color = group_type)) +
      geom_point(aes(size = career_length_yrs), alpha = .7) +
      geom_text_repel(aes(label = if_else(optimal_k > 7 | n_authors > 7, scholar_name, "")), 
                      color = "black", size = 3.2, fontface = "italic",
                      data = d_analysis) +
      theme_minimal() +
      scale_color_brewer(palette = "Dark2") +
      scale_y_continuous(breaks = 1:8) + 
      labs(x = "Average Number of Co-Authors", 
           y = "Number of Research Topics (Optimal K)", 
           color = "", size = "Career Length\n(Years)")
    
    summary(lm(optimal_k ~ n_authors + group_type, d_analysis))
    # I see no evidence of an effect
    
  # By Raw Variance
    # Men or Women, University or Institute
    d_analysis %>% 
      group_by(group_type, gender) %>% 
      summarise(total_var = mean(total_var)) %>% 
      ggplot(aes(group_type, total_var, group = gender, color = gender)) +
      geom_line(size = 1.5) +
      geom_quasirandom(width = .1, data = d_analysis, dodge.width=.4, alpha = .8) + 
      geom_text_repel(aes(label = scholar_name),
                      size = 3,
                      data = d_analysis) +
      theme_minimal() +
      theme(axis.text.y = element_blank()) +
      labs(x = "", y = "Total Variance of Publication Semantic Embeddings", color = "Gender")
    
    total_var_anova <- aov_ez(id = "scholar_name", dv = "total_var",
                              d_analysis, between = c("gender", "group_type"),
                              anova_table = list(es = "pes"))
    total_var_anova
    emmeans::emmeans(total_var_anova, c("gender", "group_type"))
    # Average Number of Co-Authors (check for interaction with above)
    d_analysis %>% 
      ggplot(aes(n_authors-1, total_var, color = group_type)) +
      geom_point(aes(size = career_length_yrs), alpha = .8) +
      geom_text_repel(aes(label = if_else(total_var > 2800 | n_authors > 7, scholar_name, "")), 
                      color = "black", size = 3.2, fontface = "italic",
                      data = d_analysis) +
      theme_minimal() +
      theme(axis.text.y = element_blank()) +
      scale_color_brewer(palette = "Dark2") +
      labs(x = "Average Number of Co-Authors", 
           y = "Total Variance of Publication Semantic Embeddings", 
           color = "", size = "Career Length\n(Years)")

    cor(d_analysis$n_authors, d_analysis$total_var)
    