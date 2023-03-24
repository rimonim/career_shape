library(tidyverse)
library(psych)      # for PCA
library(lubridate)
library(timetk)
library(tidymodels)
library(modeltime)  # for time series model building
library(cowplot)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load + Wrangle Data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

setwd('/Users/louisteitelbaum/Projects/quantitative_history/data')
filenames <- list.files()
d <- lapply(filenames, read_csv)

# remove_almog <- function(x) rename_with(x, ~str_remove(.x, 'almog_'), starts_with("almog_"))
# d <- lapply(d, remove_almog)

remove_first <- function(dat) 
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# One Tibble to Rule them All
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
d_agg <- bind_rows(d)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Graphical Analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

## Pretty Plot of PC1 by PC2 for each Scholar, Colored by Years Since First Publication
p1 <- d_agg %>% 
  # Remove outliers from the beginning of the career 
  # (A few of these are clearly mistakes by the Semantic Scholar Algorithm)
  left_join(d_agg %>% 
              group_by(scholar_name) %>% 
              summarise(min_date = quantile(date, probs = c(.01), na.rm = T))) %>% 
  filter(date > min_date) %>% 
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

d_agg %>% 
  # Remove outliers from the beginning of the career
  left_join(d_agg %>% 
              group_by(scholar_name) %>% 
              summarise(min_date = quantile(date, probs = c(.01), na.rm = T))) %>% 
  filter(date > min_date) %>% 
  # Construct Graph
  group_by(scholar_name) %>% 
  plot_time_series(date, pc1, 
                   .facet_ncol  = 4, .facet_nrow  = 5,
                   .facet_scales = "fixed",
                   .smooth = T,
                   .interactive = FALSE)

almog_simchon_pca 
yoav_kessler_pca
steven_pinker_pca
`daniel_l schacter_pca`

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Modeling
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

glmnet_mod <- linear_reg(penalty = 0.01) %>% 
  set_engine("glmnet") %>% 
  fit(
    pc1 ~ month(date, label = TRUE)
          + as.numeric(date)
          + n_citations
          + n_authors,
    training(d_agg %>% 
               filter(scholar_name == "Michael Anderson") %>% 
               time_series_split())
  )

modeltime_table(glmnet_mod)  %>% 
  modeltime_forecast(new_data = d_agg %>% 
                       filter(scholar_name == "Michael Anderson"),
                     actual_data = d_agg %>% 
                       filter(scholar_name == "Michael Anderson")
                     ) %>% 
  plot_modeltime_forecast()




  