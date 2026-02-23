####
# Turkana Miocene to modern East Africa traits against environmental proxies
# Jeremias Gloeggler 2025/5/8
####

pacman::p_load(
  "tidyverse", "ggthemes", "readxl", "palaeoverse", "ggpmisc", "ggpubr", 
  "corrplot", "deeptime"
)

destination_path <- "images/"

default_crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

theme_set(ggpubr::theme_pubclean()+theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')))

color_scheme <- c("East" = "#1C9E77",
                  "North" = "#7570B3",
                  "West" = "#D95F02",
                  "Tugen Hills"="#E72A8A",
                  "Ileret" = "#6A3D9A",
                  "Karari Ridge" = "#1F78B4",
                  "Koobi Fora Ridge"="#33A02B",
                  "Turkana - general" = "#00FFFF"
)

all_traits <- c("logBM", "HYP", "LOP", "AL", "ALX", "OL", "OLX", "GT", "BUN", "SF", "OT")

genus_list <- read_csv("genus_list.csv")

genus_list <- genus_list |> 
  filter(
# filter by redundant taxa
         is.na(Delete),
# filter to not include Homo spp.
         Genus != "Homo"|is.na(Genus)) |> #,
  mutate(Member = case_when(str_detect(Member, "Lomekwi") & !Member %in% c("Lower Lomekwi", "Upper Lomekwi") & min_ma > 3 ~ "Lower Lomekwi",
                            str_detect(Member, "Lomekwi") & !Member %in% c("Lower Lomekwi", "Upper Lomekwi") & max_ma < 3 ~ "Upper Lomekwi",
                            TRUE ~ Member
                            ))

genus_traits <- read_xlsx("genus_traits.xlsx")

# Genus trait correlation matrix ------------------------------------------

row.names(genus_traits) <- genus_traits$Genus 
genus_cor <- cor(na.omit(genus_traits[all_traits])[all_traits])

testRes.1 = cor.mtest(genus_cor, conf.level = 0.95)

corrplot(genus_cor, method = 'square',
         tl.col = 'black', cl.ratio = 0.2, 
         number.cex = 1,
         cl.cex = 1, 
         # tl.cex = .3,
         addCoef.col = 'black', tl.srt = 45,
         col = corrplot::COL2('BrBG', 100),
         # p.mat = testRes.1$p,
         insig = "blank")

# Create time bins --------------------------------------------------------

occdf <- genus_list |>
  select(Locality, Place, Member, Formation,  max_ma, min_ma) |>
  filter(!is.na(min_ma) & !is.na(max_ma)) |> 
  mutate(max_ma = if_else(max_ma-min_ma == 0, max_ma + 1e-7, max_ma)) # add 1/10 years in case of point dates, so that bin_time works with method = "majority" 

bins <- data.frame(bin = c(1:10),
                   min_ma = c(0:9),
                   max_ma = c(1:10)) |> 
  rbind(c(0, 0, 0))

binned_Mb <- bin_time(occdf = occdf,
                      bins = bins,
                      method = "majority") |>  #  "majority" gives percentage of overlap, can be used to filter out overlaps < 30-40%
  mutate(bin_assignment = case_when(Member == "modern" ~ 0, 
                                    TRUE ~ as.numeric(bin_assignment)),
         bin_midpoint = case_when(Member == "modern" ~ 0, 
                                  TRUE ~ as.numeric(bin_midpoint))
  ) |> 
  cbind(occdf$Locality) |> 
  filter(overlap_percentage > 40)

# combine species-locs-mb with traits -------------------------------------

min_body_mass_kg <- 1
min_genus_number_per_loc <- 5

loc_mb_mean_scores <- left_join(genus_list, genus_traits,
                                by = join_by(Genus)) |> 
# filter by minimum body mass
  filter(Body_Mass_Kg >= min_body_mass_kg|is.na(Genus)) |> 
# impute traits at higher taxonomic levels
  mutate(grouping_variable = coalesce(Tribe, Subfamily, Family)) |> 
  group_by(grouping_variable) |> 
  mutate(across(all_of(all_traits), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) |> 
# calculate locality means
  group_by(Locality, Place, Member) |>
  mutate(across(c(all_of(all_traits), "Body_Mass_Kg"), function(.x) mean(.x, na.rm = TRUE),
                .names = "{col}_mean"),
         across(c(all_of(all_traits), "Body_Mass_Kg"), function(.x) sd(.x, na.rm = TRUE),
                .names = "{col}_SD")) |>
# filter fossil localities by minimum number of genera
  group_by(Locality, Place, Member) |> 
  mutate(n_genus = n()) |> 
  ungroup() |> 
  filter(n_genus >= min_genus_number_per_loc|Member == "modern") |> 
  select(Region, Subregion, Locality, Place, Member, Formation, ends_with("_mean"), ends_with("_SD"), n_genus,
         min_ma, max_ma) |> 
  distinct(Locality, Place, Member, .keep_all = TRUE) |>
  inner_join(binned_Mb[c("Locality", "Place","Formation", "Member", "n_bins",
                        "bin_assignment", "bin_midpoint", "max_ma", "min_ma")],
            by = join_by(Locality, Place, Member, Formation),
            multiple = "first", suffix = c("", "_bin")) |>
  # left_join(bins, by = join_by(bin_assignment == bin), suffix = c("", "_bin")) |>
  rename_with(.fn = \(x)sub("_mean$", "", x))|> 
  filter(!is.na(min_ma)) |>  
  mutate(mean_ma = (max_ma + min_ma)/2,
         Member = str_replace_all(Member, c("TuluBor" = "Tulu Bor",
                                            "UpperBurgi" = "Upper Burgi",
                                            "lower Kalochoro" = "Kalochoro")))

# Plot temporal trends of functional traits -------------------------------

plotlist <- list()

for (trait in all_traits){
  
  p <- loc_mb_mean_scores |> 
    pivot_longer(cols = all_of(all_traits), names_to = "Trait", values_to = "Score") |> #ends_with("_SD), all_traits, all_traits
    filter(Trait == trait) |> 
    ggplot(aes(mean_ma, Score, color = Region, fill = Region))+
    geom_jitter(aes(group = paste0(mean_ma, Region)), alpha = .5)+
    stat_smooth(span = .9, se = FALSE, alpha = .8, linewidth = 1)+
    scale_color_manual(values = color_scheme, labels = c("East" = "East Turkana", 
                                                         "North" = "Lower Omo Valley", 
                                                         "West" = "West Turkana"))+
    scale_fill_manual(values = color_scheme, labels = c("East" = "East Turkana", 
                                                        "North" = "Lower Omo Valley", 
                                                        "West" = "West Turkana"))+
    scale_x_reverse(breaks = stages$min_age, labels = round(stages$min_age, 1))+
    scale_x_reverse(breaks = stages$min_age, labels = function(b) ifelse(b < 0.1, "", round(stages$min_age, 1)))+
    coord_geo(
      pos = as.list(rep("bottom", 2)),
      dat = list("stages", "epochs"),
      alpha = .5,
      height = list(unit(1.2, "lines"), unit(1.2, "lines")),
      rot = list(0, 0), size = list(6, 6), abbrv = list(TRUE, FALSE),
      center_end_labels = TRUE,
      skip = c("Holocene", "Ionian", "Greenlandian", "Northgrippian", "Meghalayan",
               "Late Pleistocene")
    )+
    labs(x = "Age (Ma)", y = trait)+
    # annotate("text", x = 7.5, y = max(loc_mb_mean_scores[[trait]])*.9, label = paste0(trait),
    #          size = 6)+
    facet_grid(rows = vars(Trait), 
               axes = "all_y",
               scale = "free")+
    theme(legend.position = "top",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          strip.background = element_blank(),
          strip.text.y = element_blank(),
          aspect.ratio = 1/3)
    
    
  p
  
  # ggsave(paste0(destination_path, trait, "_temporal_trend.svg"),
  # p, height = 6, width = 12)
  
  plotlist[[trait]] <- p
  
}

library(gridExtra)

pdf(paste0(destination_path, "all_traits_temporal_trends.pdf"), onefile = TRUE)
for (i in seq(length(plotlist))) {
  do.call("grid.arrange", plotlist[i])  
}
dev.off()

p1 <- loc_mb_mean_scores |> 
  pivot_longer(cols = all_of(all_traits), names_to = "Trait", values_to = "Score") |> #ends_with("_SD), all_traits, all_traits
  mutate(Region = factor(Region, levels = c("East", "West", "North", "Tugen Hills")),
         Trait = fct_relevel(Trait, all_traits)) |> 
  ggplot(aes(mean_ma, Score, color = Region, fill = Region))+
  stat_smooth(span = .9, se = FALSE, alpha = .8, linewidth = 1)+
  scale_color_manual(values = color_scheme, labels = c("East" = "East Turkana", 
                                                       "North" = "Lower Omo Valley", 
                                                       "West" = "West Turkana"))+
  scale_fill_manual(values = color_scheme, labels = c("East" = "East Turkana", 
                                                      "North" = "Lower Omo Valley", 
                                                      "West" = "West Turkana"))+
  scale_x_reverse(breaks = stages$min_age, labels = function(b) ifelse(b < 0.1, "", round(stages$min_age, 1)))+
  scale_y_continuous(expand = c(0,0), position = "right") +
  coord_geo(
    pos = as.list(rep("bottom", 2)),
    dat = list("stages", "epochs"),
    alpha = .5,
    height = list(unit(1.6, "lines"), unit(1.6, "lines")),
    rot = list(0, 0), size = list(8, 8), abbrv = list(TRUE, FALSE),
    center_end_labels = TRUE,
    skip = c("Holocene", "Ionian", "Greenlandian", "Northgrippian", "Meghalayan",
             "Late Pleistocene")
  )+
  labs(x = NULL, y = NULL)+
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.border = element_rect(colour = "grey40", fill=NA, linewidth=1),
        legend.text=element_text(size = 20),
        legend.key.spacing.x = unit(1, "cm"),
        legend.title = element_blank(),
        panel.spacing = unit(1, 'points'))+
  guides(colour = guide_legend(override.aes = list(linewidth = 8)))+
  facet_grid(rows = vars(Trait), scales = "free_y", switch = "y")+
  ggh4x::facetted_pos_scales(y = list(
    Trait == "logBM" ~ scale_y_continuous(breaks = c(4.5, 5.5, 6.5), position = "right"),
    Trait == "LOP" ~ scale_y_continuous(breaks = c(.95, 1.15, 1.35), position = "right"),
    Trait == "AL" ~ scale_y_continuous(breaks = c(.25, .4, .55), position = "right"),
    Trait == "GT" ~ scale_y_continuous(breaks = c(.15, .25, .35), position = "right"),
    Trait == "OT" ~ scale_y_continuous(breaks = c(.1, .25, .4), position = "right")
    ))

ggsave(paste0(destination_path, "all_traits_temporal_trend.svg"),
p1, height = 20, width = 12)


# read in environmental proxies -------------------------------------------

mammal_enamel_isotopes <- read_xlsx("large_mammal_isotopes.xlsx") |> 
  mutate(Member = case_when(Formation == "Shungura" ~ paste0(Formation, " ", str_extract(Member, "^[A-Z]")), 
                            Region == "Tugen Hills" ~ Formation,
                            str_detect(Member, "Lomekwi") & !Member %in% c("Lower Lomekwi", "Upper Lomekwi") & mean_ma > 3 ~ "Lower Lomekwi",
                            str_detect(Member, "Lomekwi") & !Member %in% c("Lower Lomekwi", "Upper Lomekwi") & mean_ma < 3 ~ "Upper Lomekwi",
                            TRUE ~ Member),
         Member = str_replace(Member, "middle", "Middle"),
         bin_assignment = ifelse(mean_ma == 0, 0, cut(mean_ma, breaks = seq(0, 10),
                                                      labels = seq(1,10)))) |> 
  group_by(Region, Member) |> 
  filter(n() >= 5) |> 
  ungroup()

pedogenic_isotopes <-  read_xlsx("pedogenic_carbonate_isotopes.xlsx") |> 
  mutate(Member = case_when(
                            Formation == "Tugen Hills - modern" ~ "modern",
                            Region == "Tugen Hills" ~ Formation,
                            str_detect(Member, "Lomekwi") & !Member %in% c("Lower Lomekwi", "Upper Lomekwi") & mean_ma > 3 ~ "Lower Lomekwi",
                            str_detect(Member, "Lomekwi") & !Member %in% c("Lower Lomekwi", "Upper Lomekwi") & mean_ma < 3 ~ "Upper Lomekwi",
                            TRUE ~ Member),
         Member = str_replace(Member, "middle", "Middle"),
         Formation = if_else(Formation == "Tugen Hills - modern", "modern", Formation),
         bin_assignment = ifelse(mean_ma == 0, 0, cut(mean_ma, breaks = seq(0, 10),
                                                      labels = seq(1,10))))|> 
  group_by(Region, Member) |> 
  filter(n() >= 5) |> 
  ungroup()

region_mb_mean_enamel_isotope <- mammal_enamel_isotopes |>
  filter(!(Member == "modern" & Region == "West")) |>
  mutate(Member = case_when(Formation == "Kanapoi" ~ "Kanapoi", 
                            TRUE ~ Member),
         mean_d13C = mean(d13C_VPDB, na.rm = TRUE),
            SD_d13C = sd(d13C_VPDB, na.rm = TRUE),
            mean_d18O = mean(d18O_VPDB, na.rm = TRUE),
            SD_d18O = sd(d18O_VPDB, na.rm = TRUE),
            .by = c(Region, Member)) |> 
  dplyr::select(Region, Member, Formation, mean_d13C, SD_d13C, mean_d18O, SD_d18O) |> 
  distinct(Region, Member, Formation, .keep_all = TRUE)

region_mb_mean_CO_isotope <- pedogenic_isotopes |> 
  mutate(mean_d13C = mean(d13C, na.rm = TRUE),
         SD_d13C = sd(d13C, na.rm = TRUE),
         mean_d18O = mean(d18O, na.rm = TRUE),
         SD_d18O = sd(d18O, na.rm = TRUE),
         .by = c(Region, Member)) |> 
  select(Region, Member, Formation, mean_d13C, SD_d13C, mean_d18O, SD_d18O) |> 
  distinct(Region, Member, .keep_all = TRUE)

bin_mean_trait <- loc_mb_mean_scores |>
  group_by(Member, Formation) |> 
  filter(max(max_ma) - min(min_ma) < 2) |> 
  ungroup() |> 
  mutate(across(c(all_of(all_traits), "Body_Mass_Kg"), function(.x) mean(.x, na.rm = TRUE),
                .names = "mean_{col}"),
         across(c(all_of(all_traits), "Body_Mass_Kg"), function(.x) sd(.x, na.rm = TRUE),
                .names = "SD_{col}"),
         .by = c(bin_assignment))|> 
  select(Region, Member, Formation, bin_assignment, starts_with("mean_"), starts_with("SD_"), mean_ma, min_ma, max_ma) |> 
  distinct(bin_assignment, .keep_all = TRUE)

region_mb_mean_trait <- loc_mb_mean_scores |> 
  group_by(Member) |> 
  filter(max(max_ma) - min(min_ma) < 2) |>
  ungroup() |> 
  mutate(across(c(all_of(all_traits), "Body_Mass_Kg"), function(.x) mean(.x, na.rm = TRUE),
                   .names = "mean_{col}"),
            across(c(all_of(all_traits), "Body_Mass_Kg"), function(.x) sd(.x, na.rm = TRUE),
                   .names = "SD_{col}"),
            .by = c(Region, Member))|> 
  select(Region, Subregion, Member, Formation, starts_with("mean_"), starts_with("SD_"), mean_ma, min_ma, max_ma,
         bin_assignment) |> 
  distinct(Member, Region, .keep_all = TRUE)


# modern trait scenario ---------------------------------------------------

region_mb_mean_trait.modern <- loc_mb_mean_scores |> 
  group_by(Member) |> 
  filter(max(max_ma) - min(min_ma) < 2,
         Member == "modern" & str_detect(Locality, "HID") & Region == "Tugen Hills"| # Sibiloi HID: 56453
         Member != "modern") |>
  ungroup() |> 
  mutate(across(c(all_of(all_traits), "Body_Mass_Kg"), function(.x) mean(.x, na.rm = TRUE),
                .names = "mean_{col}"),
         across(c(all_of(all_traits), "Body_Mass_Kg"), function(.x) sd(.x, na.rm = TRUE),
                .names = "SD_{col}"),
         .by = c(Region, Member))|> 
  select(Locality, Region, Subregion, Member, Formation, starts_with("mean_"), starts_with("SD_"), mean_ma, min_ma, max_ma,
         bin_assignment) |> 
  group_by(Region, Member) |>
  filter(row_number() == 1) |> 
  ungroup()

region_mb_combined_enamel <- inner_join(region_mb_mean_trait,
                                 region_mb_mean_enamel_isotope,
                                 by = join_by(Region, Member, Formation),
                                 na_matches = "never")

region_mb_combined_enamel.modern <- inner_join(region_mb_mean_trait.modern,
                                         region_mb_mean_enamel_isotope,
                                         by = join_by(Region, Member, Formation),
                                         na_matches = "never")

region_mb_combined_CO <- inner_join(region_mb_mean_trait,
                                    region_mb_mean_CO_isotope,
                                    by = join_by(Region, Member, Formation),
                                    na_matches = "never")

region_mb_combined_CO.modern <- inner_join(region_mb_mean_trait.modern,
                                     region_mb_mean_CO_isotope,
                                     by = join_by(Region, Member, Formation),
                                     na_matches = "never")

region_bin_combined_total <- inner_join(region_mb_combined_CO, region_mb_combined_enamel[c("Region", "Member", "mean_d13C")],
                                        by = join_by(Region, Member),
                                        suffix = c("_CO", "_enamel"))

region_bin_combined_total.modern <- inner_join(region_mb_combined_CO.modern, region_mb_combined_enamel.modern,
                                        by = join_by(Region, Member),
                                        suffix = c("_CO", "_enamel"))
### 
# region_mb_combined_CO |>
#   # filter(Member != "modern") |>
#   GGally::ggscatmat(columns = paste0("mean_", all_of(c(all_traits, c("d13C", "d18O")))) , color = "Region", corMethod = "spearman") +
#   scale_color_manual(values = color_scheme)+
#   ggpubr::theme_pubclean()+
#   labs(x = "X", y = "Y")+
#   theme(strip.background = element_rect(fill = "grey90"),
#         axis.text = element_text(size = 8),
#         legend.position = "right",
#         legend.key = element_blank())

# Scatter plots with correlations -----------------------------------------

corr_isotope <- "d18O"
isotope_labels <- c("d13C" = "*&delta;*<sup>13</sup>C<br />(&permil; VPDB)", 
                    "d18O" = "*&delta;*<sup>18</sup>O<br />(&permil; VPDB)",
                    "woody_cover" = "Woody cover (%)")[corr_isotope]

time_bin_cutoff <- 2.7

traits_cor_dataset <- region_mb_combined_CO.modern # region_mb_combined_enamel or region_mb_combined_CO

dataset_name_string <- ifelse(identical(traits_cor_dataset, region_mb_combined_enamel), "enamel", "CO")

for (trait in all_traits[2]){
  mean_correlation_plot_binned <- traits_cor_dataset |> 
    mutate(time_bin = mean_ma > time_bin_cutoff,
           time_bin = if_else(time_bin, paste0("> ", time_bin_cutoff," Ma"), 
                              paste0("<= ", time_bin_cutoff," Ma"))) |> 
    # filter(Member != "modern") |>
    ggplot(aes_string(paste0("mean_", corr_isotope), paste0("mean_", trait)))+
    geom_point(aes(color = Region, size = as.factor(bin_assignment)))+
    geom_errorbarh(aes_string(xmin = paste0("mean_", corr_isotope,"-SD_", corr_isotope), 
                              xmax = paste0("mean_", corr_isotope,"+SD_", corr_isotope), 
                              y = paste0("mean_", trait), color = "Region",
                              height = 0.01),
                   alpha = .3)+
    geom_errorbar(aes_string(ymin = paste0("mean_", trait, "- SD_", trait), 
                             ymax = paste0("mean_", trait, "+ SD_", trait),
                             x = paste0("mean_", corr_isotope), color = "Region"),
                  alpha = .3)+   
    stat_cor(color = "grey30", label.y.npc = "bottom")+
    geom_smooth(color = "grey30", method = "lm", se = FALSE, linetype = "dashed")+
    scale_color_manual(values = color_scheme)+
    ggnewscale::new_scale_color()+
    stat_cor(aes(color = time_bin))+ # labels: https://cran.r-project.org/web/packages/ggpmisc/vignettes/model-based-annotations.html
    geom_smooth(aes(color = time_bin), method = "lm", se = FALSE)+
    scale_color_viridis_d("Time bin", end = .5)+
    labs(color = "Region", size = "Time bin",
         y = trait)+
    xlab(isotope_labels)+
    guides(color = guide_legend(override.aes = list(size = 5)),
           size = guide_legend(ncol=2))+
    theme(axis.title.x = ggtext::element_markdown(),
          legend.position = "right")
  
# ggsave(paste0(destination_path, trait, "_trait_", dataset_name_string, "_", corr_isotope,"_plot_", time_bin_cutoff, "Ma.svg"), width = 8, height = 6)
  
}


plotly::ggplotly(mean_correlation_plot_binned)
# Correlation matrices ----------------------------------------------------

time_bin_cutoff.1 <- 5.3
time_bin_cutoff.2 <- 2.7

trait_metrics <- "mean"
isotope_metrics <- "mean"

cor_dataframe.1 <- region_mb_combined_CO |>
  # filter(mean_ma <= 2.7 | mean_ma >= 4) |>
  drop_na(paste0(c("mean" = "mean_", "SD" = "SD_")[trait_metrics], all_traits), 
          paste0(c("mean" = "mean_", "SD" = "SD_")[isotope_metrics], c("d13C", "d18O"))) |>  #   drop_na(paste0("SD_", all_traits), "SD_d13C", "SD_d18O")
  filter(mean_ma <= time_bin_cutoff.2, Member != "modern")

cor_dataframe.2 <-  region_mb_combined_enamel |>
  # filter(mean_ma <= 2.7 | mean_ma >= 4) |>
  drop_na(paste0(c("mean" = "mean_", "SD" = "SD_")[trait_metrics], all_traits), 
          paste0(c("mean" = "mean_", "SD" = "SD_")[isotope_metrics], c("d13C"))) |>  #   drop_na(paste0("SD_", all_traits), "SD_d13C", "SD_d18O")
  filter(mean_ma <= time_bin_cutoff.2, Member != "modern")

traits.1 <- cor_dataframe.1[paste0(c("mean" = "mean_", "SD" = "SD_")[trait_metrics], all_traits)]
isotopes.1 <- cor_dataframe.1[paste0(c("mean" = "mean_", "SD" = "SD_")[isotope_metrics], c("d13C", "d18O"))]

traits.2 <- cor_dataframe.2[paste0(c("mean" = "mean_", "SD" = "SD_")[trait_metrics], all_traits)]
isotopes.2 <- cor_dataframe.2[paste0(c("mean" = "mean_", "SD" = "SD_")[isotope_metrics], c("d13C"))]

cor_matrix.1 <- cor(cbind(traits.1, isotopes.1), method = "pearson") # pearson pr spearman
cor_matrix.2 <- cor(cbind(traits.2, isotopes.2), method = "pearson") # pearson pr spearman

testRes.1 = cor.mtest(cor_dataframe.1[c(paste0(c("mean" = "mean_", "SD" = "SD_")[trait_metrics], all_traits), 
                                        paste0(c("mean" = "mean_", "SD" = "SD_")[isotope_metrics], c("d13C", "d18O")))], conf.level = 0.95)
testRes.2 = cor.mtest(cor_dataframe.2[c(paste0(c("mean" = "mean_", "SD" = "SD_")[trait_metrics], all_traits), 
                                        paste0(c("mean" = "mean_", "SD" = "SD_")[isotope_metrics], c("d13C")))], conf.level = 0.95)

par(mar = rep(3, 4))
par(oma = rep(10, 4))

corrplot(cor_matrix.1*100, method = 'square', type = 'upper',
         tl.col = 'black', cl.ratio = 0.2, 
         number.cex = .9,
         number.digits = 0,
         cl.cex = 1, 
         cl.pos = FALSE,
         # tl.cex = .3,
         addCoef.col = 'black', tl.srt = 45, tl.pos = "t",
         col = corrplot::COL2('BrBG', 100), 
         # p.mat = testRes.1$p,
         is.corr = FALSE,
         # diag = FALSE,
         insig = "blank")

corrplot(cor_matrix.2*100, method = 'square', type = 'lower', 
         tl.col = 'black', cl.ratio = 0.2, 
         number.cex = .9,
         number.digits = 0,
         cl.cex = 1,
         cl.pos = FALSE,
         # tl.cex = .3,
         addCoef.col = 'black', tl.srt = 0, 
         col = corrplot::COL2('BrBG', 100), tl.pos = "l",
         # p.mat = testRes.2$p,
         insig = "blank",
         is.corr = FALSE,
         diag = F,
         add = TRUE)

# height = 1001, width = 1400

# Canonical correlation analysis (CCorA) ----------------------------------
ccora_dataframe <- region_mb_combined_CO |>                 # region_mb_combined_CO and region_mb_combined_enamel
    filter(Region != "Tugen Hills")
trait_metrics <- "mean"
isotope_metrics <- "mean"

traits <- ccora_dataframe[paste0(c("mean" = "mean_", "SD" = "SD_")[trait_metrics], all_traits)]
isotopes <- ccora_dataframe[paste0(c("mean" = "mean_", "SD" = "SD_")[isotope_metrics], c("d13C", "d18O"))]

cc <- yacca::cca(traits, isotopes, na.rm = TRUE)

par(mar = c(5.1, 4.1, 4.1, 2.1))
par(oma = rep(0, 4))
yacca::helio.plot(cc, # save with ratio 600:660
                  xvlab = str_replace_all(cc$xlab, c("mean_" = "", "_CO" = "")),
                  x.name = "Functional\ntraits",
                  yvlab = str_replace_all(cc$ylab, c("mean_" = "", "d" = paste0("\u03B4"), "13" = "¹³",
                                                     "18" = "¹⁸")),
                  y.name = "Environmental\nproxies",
                  main = NULL, 
                  wid.fac = 1,
                  type = "correlation" # correlation or variance
                  )

plot(candisc::cancor(traits, isotopes), #ellipse.args = list(levels = 0.95),
     labels = str_replace(ccora_dataframe$Member, " ", " \n"),
     pch = 16,
     id.n = 2)

# Test all possible trait combinations ------------------------------------

xcombs <- unlist(
  lapply(1:length(traits), function(k) combn(names(traits), k, simplify = FALSE)),
  recursive = FALSE
)

cca_results <- lapply(xcombs, function(comb) {
  X <- as.matrix(ccora_dataframe[, comb, drop = FALSE])
  Y <- as.matrix(ccora_dataframe[, names(isotopes), drop = FALSE])
  fit <- yacca::cca(X, Y)
  list(
    fit = fit,
    cancor = fit$corr,
    RI = fit$yvrd
  )
})

names(cca_results) <- sapply(xcombs, paste, collapse = "_")

summary_table <- data.frame(
  model = names(cca_results),
  cancor1 = sapply(cca_results, function(x) x$cancor[1]),
  RI1 = sapply(cca_results, function(x) x$RI[1])
) |> 
  arrange(desc(RI1), desc(cancor1))

CCA::plt.cc(CCA::cc(traits, isotopes), type = "b", var.label = TRUE,
            ind.names = paste(str_replace(ccora_dataframe$Member, " ", " \n"), sep = "\n"))

ggplot()+
  geom_line(aes(x = 0, y = c(-Inf, Inf)))+
  geom_line(aes(y = 0, x = c(-Inf, Inf)))+
  ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = 1), color = "black")+
  ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = .5), color = "black", linetype = "dashed")+
  geom_text(aes(cc[["xstructcorr"]][,1], cc[["xstructcorr"]][,2], label = paste("mean_", all_traits)))+
  geom_segment(aes(x = 0, y = 0, xend = cc[["xstructcorr"]][,1], yend = cc[["xstructcorr"]][,2]), arrow = arrow(length = unit(0.5, 'cm')))+
  geom_text(aes(cc[["ycrosscorr"]][,1], cc[["ycrosscorr"]][,2], label = c("mean_d13C", "mean_d18O")),
            color = "red")+
  geom_segment(aes(x = 0, y = 0, xend = cc[["ycrosscorr"]][,1], yend = cc[["ycrosscorr"]][,2]), arrow = arrow(length = unit(0.5, 'cm')),
               color = "red")+
  ggrepel::geom_text_repel(aes(cc[["canvarx"]][,1]/max(cc[["canvarx"]][,1]), cc[["canvarx"]][,2]/max(cc[["canvarx"]][,2]), label = str_replace(ccora_dataframe$Member, " ", " \n"),
                color = ccora_dataframe$Region), bg.color = "white", nudge_x = 0, nudge_y = 0)+
  scale_color_manual(values = color_scheme, "Region")+
  labs(x = "Dim 1", y = "Dim 2")+
  coord_fixed()+
  theme_test()

# check significance of canonical variates
CCP::p.asym(cc$corr, dim(traits)[1], length(traits), length(isotopes), tstat = "Wilks") # corresponds to yacca::F.test.cca
CCP::p.asym(cc$corr, dim(traits)[1], length(traits), length(isotopes), tstat = "Hotelling")
CCP::p.asym(cc$corr, dim(traits)[1], length(traits), length(isotopes), tstat = "Pillai")
CCP::p.asym(cc$corr, dim(traits)[1], length(traits), length(isotopes), tstat = "Roy")

# only first CV is significant

CV1_X <- as.matrix(traits) %*% cc$xcoef[, 1]
CV1_Y <- as.matrix(isotopes) %*% cc$ycoef[, 1]

cca_df <- ccora_dataframe |> 
  mutate(CV1_X = CV1_X,
         CV1_Y = CV1_Y,
         anomaly = CV1_Y/CV1_X - mean(CV1_Y/CV1_X))

# Plot canonical correlation analysis results -----------------------------

cca_df$bin_assignment <- as.factor(cca_df$bin_assignment)

p1 <- cca_df |> 
  ggplot(aes(x = bin_assignment, y = CV1_X,  group = bin_assignment))+
  geom_boxplot(width=0.6, outlier.shape = NA, color = "grey30", fill = "grey30", alpha = .2)+
  ggnewscale::new_scale_color()+
  geom_jitter(width=0.15, size = 2, aes(color = Region))+
  scale_color_manual(values = color_scheme)+
  theme(legend.position="none")+
  scale_x_discrete(limits = rev, "Time bin")+
  theme(panel.border = element_rect(colour = "grey30", fill=NA, linewidth = .5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15))

p2 <- cca_df |> 
  ggplot(aes(x = bin_assignment, y = CV1_Y, group = bin_assignment))+
  geom_boxplot(width=0.6, outlier.shape = NA, color = "grey30", fill = "grey30", alpha = .2)+
  ggnewscale::new_scale_color()+
  geom_jitter(width=0.15, size = 2, aes(color = Region))+
  scale_color_manual(values = color_scheme)+
  theme(legend.position = "none")+
  scale_x_discrete(limits = rev, "Time bin")+
  scale_y_continuous(position = "right")+
  theme(panel.border = element_rect(colour = "grey30", fill=NA, linewidth = .5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15))

ggpubr::ggarrange(p1, p2) # 500:500

main.plot <- cca_df |> # height: 700 - width: 550
  # filter(mean_ma <= 2) |> 
  ggplot(aes(x = CV1_X,y = CV1_Y, color = Region, size = Member))+
  geom_line(stat="smooth", method = "lm",
            linewidth = 1,
            linetype ="dashed",
            alpha = 0.5)+
  geom_point()+
  ggrepel::geom_text_repel(aes(label = str_replace_all(Member, c(" " = " \n",
                                                                 "1" = "I",
                                                                 "2" = "II",
                                                                 "3" = "III",
                                                                 "4" = "IV"))), 
                           
                           bg.color = "white")+ # text: bin_assignment or str_replace(Member, " ", " \n")
  ggfocus::scale_size_focus("modern", size_focus = 5, size_other = 4)+
  # stat_cor()+
  # stat_poly_eq(use_label(c("eq", "R2")))+
  scale_y_continuous(limits = c(min(CV1_Y-0.05), max(CV1_Y+0.05)))+
  scale_color_manual(values = color_scheme)+
  guides(colour = guide_legend(override.aes = list(size = 5)),
         )+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15), 
        legend.position = "top")

main.plot

cca_anomaly_scatter <- cca_df |> # 550:700
  ggplot(aes(mean_ma, anomaly))+
  geom_ribbon(color = "grey",
              stat = "smooth",
              se = TRUE,
              alpha = 0, # or, use fill = NA
              linetype = "dotted")+
  geom_smooth(se = FALSE, color = "grey", span = .9)+
  geom_segment(aes(xend = mean_ma, yend = 0), alpha = .2)+
  geom_point(aes(color = Region), size = 4)+
  geom_abline(slope = 0, alpha = .2)+
  scale_color_manual(values = color_scheme)+
  scale_x_reverse()+
  # scale_y_continuous(breaks = seq(-.5,.3,.2))+
  labs(x = "Age (Ma)", y = "Scaled CV ratio")+
  theme_classic()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15), 
        legend.position = "top")

cca_anomaly_scatter

shapiro.test(cca_df$anomaly)

outliers::grubbs.test(cca_df$anomaly)
outliers::grubbs.test(cca_df$anomaly, type = 20)
outliers::grubbs.test(cca_df$anomaly, opposite = TRUE)
EnvStats::rosnerTest(cca_df$anomaly, k = 5)

q3 <- quantile(cca_df$anomaly, 0.75)
iqr <- IQR(cca_df$anomaly)
upper_bound <- q3 + 1.5*iqr

outliers_upper <- cca_df$anomaly[cca_df$anomaly > upper_bound]
outliers_upper

q1 <- quantile(cca_df$anomaly, 0.25)
iqr <- IQR(cca_df$anomaly)
lower_bound <- q1 - 1.5*iqr

outliers_lower <- cca_df$anomaly[cca_df$anomaly < lower_bound]
outliers_lower

anomaly_boxplot <- cca_df |> 
  ggplot(aes(y = anomaly)) + 
  geom_boxplot(color = "grey")

anomaly_boxplot

combined_anomaly_plot <- cca_anomaly_scatter |> 
  aplot::insert_right(anomaly_boxplot + theme_void(), width = .2)
combined_anomaly_plot # 550:700
