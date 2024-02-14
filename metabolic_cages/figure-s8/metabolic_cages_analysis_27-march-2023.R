
# Import libraries --------------------------------------------------------

library(pzfx)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(svglite)

## Define Functions ---------------------------------------------------------------

### Data Transformations -----------------------------------------
# A function to replace commas for decimal points and convert to number
sub_and_num <- function(x){
  x_pt <- str_replace(x, c(","), ".")
  x_num <- as.numeric(x_pt)
}

# A function to extract table from table list
genotype <- function(name){
  if (str_detect(name, 'WT')){
    return('WT')
  }else{
    return('KO')
  }
}

convert_table <- function(table_list, table_num){
  df <- table_list[[table_num]] %>% mutate_at(2:11, sub_and_num)
  df_t <- as_tibble(t(df[, -1]))
  names(df_t) <- df[, 1]
  rownames(df_t) <- names(df)[-1]
  df_t['genotype'] <- unlist(lapply(rownames(df_t), genotype))
#  df_t['genotype'] <- factor(df_t['genotype'], levels = c("WT", "KO"))
  return(df_t)
}

# A function that properly merges cd and hfd tables
merge_cd_hfd <- function(cd_tables, hfd_tables, table_num){
  cd_df <- convert_table(cd_tables, table_num)
  cd_df['diet'] <- rep('CD', nrow(cd_df))
  hfd_df <- convert_table(hfd_tables, table_num)
  hfd_df['diet'] <- rep('HFD', nrow(hfd_df))
  df_full <- rbind(cd_df, hfd_df)
  print(rownames(df_full))
  names_hfd <- paste0('h', row.names(hfd_df))
  rnames_full <- c(row.names(cd_df), names_hfd)
  rownames(df_full) <- rnames_full
  return(df_full)
}

# A function to transform data for bar/dot phased plots
transform_for_barplot <- function(df, parameter,
                                  diet_subset, phase_or_composition){
  # 1 - phase, 2 - composition
  df$genotype <- factor(df$genotype, levels = c("WT", "KO"))
  phase_labels <- c("Dark Phase", "Light Phase", "Mean")
  composition_labels <- c("Fat Mass", "Fluid", "Lean Mass")
  df_labels <- data.frame(phase_labels,
                          composition_labels)
  
  df_diet_subset <- df[which(df$diet == diet_subset), ]
  
  df_light_diet_subset <- df_diet_subset[,c(2,4,5)]
  names(df_light_diet_subset)[1] <- parameter
  df_light_diet_subset$phase <- rep(df_labels[2, phase_or_composition],
                                    nrow(df_light_diet_subset))
  
  df_dark_diet_subset <- df_diet_subset[,c(1,4,5)]
  names(df_dark_diet_subset)[1] <- parameter
  df_dark_diet_subset$phase <- rep(df_labels[1, phase_or_composition],
                                   nrow(df_dark_diet_subset))
  
  df_mean_diet_subset <- df_diet_subset[,c(3,4,5)]
  names(df_mean_diet_subset)[1] <- parameter
  df_mean_diet_subset$phase <- rep(df_labels[3, phase_or_composition],
                                   nrow(df_mean_diet_subset))
  
  df_mean_diet_subset <- df_mean_diet_subset %>%
    mutate_at(1, sub_and_num)
  df_diet_subset_barplot <- rbind(df_dark_diet_subset,
                                  df_light_diet_subset, df_mean_diet_subset)
  df_diet_subset_barplot <- df_diet_subset_barplot %>%
    mutate_at(1, sub_and_num)
  return(df_diet_subset_barplot)
}

### Plotting Functions -----------------------------------------
# A function to plot phase plots
plot_phase_barplot <- function(df, df_mean, var_num,
                               scale_factor,
                               min_val = 0, max_val = 0){
  
  var_list <- list(quote(water_consumption),
                   quote(mass),
                   quote(breaks),
                   quote(kcal),
                   quote(rq))
  if (scale_factor != 0){
    overhead <- max(df_mean$mean) + range(df_mean$mean)/scale_factor
  } else {
    overhead <- max_val
  }
  ggplot(data = df_mean, aes(x = phase, y = mean,
                             group = genotype)) +
    geom_col(aes(x = phase),
             position = position_dodge(), colour = "black",
             fill = "white",
             size = 0.8) +
    geom_jitter(data = df, aes(x = phase,
                               y = eval(var_list[[var_num]]),
                               fill = genotype),
                na.rm = T,
                position = position_jitterdodge(dodge.width = 1,
                                                jitter.width = 0.8),
                shape = 21, size = 3,
                colour = "black", stroke = 1.5,
                inherit.aes = F) +
    scale_fill_manual(values = c(WT = "#B9B9B9", KO = "#009193")) +
    geom_errorbar(aes(ymin = mean-(sd/sqrt(6)),
                      ymax = mean+(sd/sqrt(6))),
                  position = position_dodge(0.9),
                  colour = "grey",
                  alpha = 1,
                  width = 0.3,
                  size = 0.5) +
    theme_classic() +
    theme(legend.position = "none") +
    scale_y_continuous(expand = c(0,0), limits = c(min_val, overhead)) +
    scale_x_discrete(expand = c(0,0))
}

# A function to plot time-series plots
plot_time_series <- function(df, diet, use_scaling = F, min_y_val = 0,
                             max_y_val = 0, scale_factor = 1){
  
  df_cd_wt <- df[which(df$genotype=='WT' & df$diet=='CD'), ]
  df_cd_ko <- df[which(df$genotype=='KO' & df$diet=='CD'), ]
  df_hfd_wt <- df[which(df$genotype=='WT' & df$diet=='HFD'), ]
  df_hfd_ko <- df[which(df$genotype=='KO' & df$diet=='HFD'), ]
  df_cd <- df[which(df$diet=='CD'), ]
  df_hfd <- df[which(df$diet=='HFD'), ]
  
  df_cd_wt <- df_cd_wt[,1:24]
  df_cd_wt <- as.data.frame(t(df_cd_wt))
  df_cd_wt <- df_cd_wt %>% mutate(mean=rowMeans(df_cd_wt))
  df_cd_wt <- df_cd_wt %>% transmute(rowwise(df_cd_wt),
                                     se=sd(c_across(1:5))/sqrt(5))
  df_cd_wt$genotype <- rep("WT", nrow(df_cd_wt))
  print(df_cd_wt)
  
  df_cd_ko <- df_cd_ko[,1:24]
  df_cd_ko <- as.data.frame(t(df_cd_ko))
  df_cd_ko <- df_cd_ko %>% mutate(mean=rowMeans(df_cd_ko))
  df_cd_ko <- df_cd_ko %>% transmute(rowwise(df_cd_ko),
                                     se=sd(c_across(1:4))/sqrt(4))
  df_cd_ko$genotype <- rep("KO", nrow(df_cd_ko))
  print(df_cd_ko)
  
  df_cd <- rbind(df_cd_wt[,6:8], df_cd_ko[,5:7])
  print(df_cd)
  
  df_hfd_wt <- df_hfd_wt[,1:24]
  df_hfd_wt <- as.data.frame(t(df_hfd_wt))
  df_hfd_wt <- df_hfd_wt %>% mutate(mean=rowMeans(df_hfd_wt))
  df_hfd_wt <- df_hfd_wt %>% transmute(rowwise(df_hfd_wt),
                                       se=sd(c_across(1:6))/sqrt(6))
  df_hfd_wt$genotype <- rep("WT", nrow(df_hfd_wt))
  print(df_hfd_wt)
  
  
  df_hfd_ko <- df_hfd_ko[,1:24]
  df_hfd_ko <- as.data.frame(t(df_hfd_ko))
  df_hfd_ko <- df_hfd_ko %>% mutate(mean=rowMeans(df_hfd_ko))
  df_hfd_ko <- df_hfd_ko %>% transmute(rowwise(df_hfd_ko),
                                       se=sd(c_across(1:8))/sqrt(8))
  df_hfd_ko$genotype <- rep("KO", nrow(df_hfd_ko))
  print(df_hfd_ko)
  
  df_hfd <- rbind(df_hfd_wt[,7:9], df_hfd_ko[,9:11])
  print(df_hfd)
  plot_data <- df_cd
  if (diet == "hfd"){plot_data <- df_hfd}
  if (min_y_val == 0 && max_y_val == 0){
    min_y_val <- min(plot_data$mean) - 10/scale_factor
    max_y_val <- max(plot_data$mean) + 10/scale_factor
  }
  
  x_time <- factor(row.names(df_cd_wt))
  x_time <- c(x_time, x_time)
  print(x_time)
  p <- ggplot(plot_data, mapping = aes(x = x_time,
                                       y = mean,
                                       group = genotype,
                                       colour = genotype,
                                       fill = genotype)) +
    geom_errorbar(aes(ymin = mean - se,
                      ymax = mean + se),
                  colour = "grey",
                  alpha = 0.6) +
    geom_line(size = 1) +
    geom_point(size = 4) +
    geom_point(shape = 1, size = 4, colour = "black", stroke = 2) +
    scale_x_discrete(limits = x_time[1:24],
                     breaks = c("19h", "22h", "1h", "4h", "7h",
                                "10h", "13h", "16h")) +
    theme_classic(base_size = 20) +
    scale_color_manual(values = c(WT = "#B9B9B9", KO = "#009193")) +
    scale_fill_manual(values = c(WT = "#B9B9B9", KO = "#009193")) +
    theme(legend.position = "none")
  if (use_scaling == T){
    p + scale_y_continuous(limits = c(min_y_val, max_y_val))
  }
}



# Preprocessing -----------------------------------------------------------

setwd("/Users/artemii/Yandex.Disk.localized/Inserm/Exps/metabolic-cages-feb-2021/metabolic-cages-2021")

PZFX_FILE_CD = '0221_RORaFXRmales_RM1_edited.pzfx'
PZFX_FILE_HFD = '0221_RORaFXRmales_HF_edited.pzfx'

tables <- pzfx_tables(PZFX_FILE_HFD)

# Old merge function
# merge_cd_hfd <- function(table_cd, table_hfd){
#   table_cd_proc <- as.data.frame(t(table_cd))
#   table_cd_proc <- table_cd_proc[-1,]
#   colnames(table_cd_proc) <- table_cd$ROWTITLE
#   table_cd_proc['diet'] <- rep('CD', nrow(table_cd_proc))
#   table_hfd_proc <- as.data.frame(t(table_hfd))
#   colnames(table_hfd_proc) <- table_hfd$ROWTITLE
#   table_hfd_proc <- table_hfd_proc[-1,]
#   table_hfd_proc['diet'] <- rep('HFD', nrow(table_hfd_proc))
#   table_merged <- rbind(table_cd_proc, table_hfd_proc)
#   return(table_merged)
# }

# Read in the data frames and apply the merge_cd_hfd() function
cd_tables <- lapply(pzfx_tables(PZFX_FILE_CD),
                    function(x) read_pzfx(path = PZFX_FILE_CD, table = x))
hfd_tables <- lapply(pzfx_tables(PZFX_FILE_HFD),
                     function(x) read_pzfx(path = PZFX_FILE_HFD, table = x))
# merged_tables <- mapply(merge_cd_hfd, cd_tables, hfd_tables, SIMPLIFY = FALSE)


# Subsetting data ---------------------------------------------------------
# > tables
# [1] "RER 24h"                       "EE 24h"                        "Food intake cum 24h"          
# [4] "XAmb 24h"                      "Z 24h"                         "RER half bars"                
# [7] "EE half bars"                  "XAmb half bars"                "Z half bars"                  
# [10] "Food intake half bars"         "Water intake half bars"        "BW"                           
# [13] "Body composition"              "Body composition (normalized)" "E Balance half bars"          
# [16] "Food intake non cum 24h"     

## Energy expenditure vs BW -------

# Merging, removing NAs
EE_df <- merge_cd_hfd(cd_tables, hfd_tables, 2)
EE_df <- na.omit(EE_df)
EE_df <- EE_df %>% mutate_at(1:24, sub_and_num)

# Summing energy expenditure
EE_df['EE_24h'] <- rowSums(EE_df[,c(-25, -26)])

# Adding body weight
BW_df <- rbind(cd_tables[[12]][1:5,], hfd_tables[[12]])
BW_col <- c(BW_df[1:5,]$WT, BW_df[1:5,]$KO, BW_df[6:13,]$WT, BW_df[6:13,]$KO)
EE_df['BW'] <- sub_and_num(BW_col[c(-8, -17, -18)])

## RER -------
RER_df <- merge_cd_hfd(cd_tables, hfd_tables, 1)
RER_df <- na.omit(RER_df)
RER_df <- RER_df %>% mutate_at(1:24, sub_and_num)

## Food intake -------
food <- merge_cd_hfd(cd_tables, hfd_tables, 3)
food <- na.omit(food)
food <- food %>% mutate_at(1:24, sub_and_num)

## Body weight -------
BW_df_num <- BW_df %>% mutate_at(1:2, sub_and_num)
bw_vector <-  c(BW_df_num[,1], BW_df_num[,2])

BW <- bw_vector
genotype_vector <- factor(c(rep("WT", 13), rep("KO", 13)),
                          levels = c("WT", "KO"))
diet <- factor(c(rep("CD", 5), rep("HFD", 8),
                 rep("CD", 5), rep("HFD", 8)),
               levels = c("CD", "HFD"))
BW_df_barplot <- data.frame(BW, genotype_vector, diet) 


mean_bw <- c(mean(BW_df_barplot$BW[1:5]),
             mean(BW_df_barplot$BW[6:13], na.rm = T),
             mean(BW_df_barplot$BW[14:18]),
             mean(BW_df_barplot$BW[19:26]))
se <- c(sd(BW_df_barplot$BW[1:5])/sqrt(5),
        sd(BW_df_barplot$BW[6:13], na.rm = T)/sqrt(6),
        sd(BW_df_barplot$BW[14:18])/sqrt(5),
        sd(BW_df_barplot$BW[19:26])/sqrt(8)) 
genotype_vector_mean <- factor(c(rep("WT", 2), rep("KO", 2)),
                               levels = c("WT", "KO"))
diet_mean <- factor(c("CD", "HFD", "CD", "HFD"),
                    levels = c("CD", "HFD"))
BW_mean_barplot <- data.frame(mean_bw, genotype_vector_mean,
                              diet_mean, se) 

## Phased (Water consumption and alike) ------

### CD ------
water_df <- merge_cd_hfd(cd_tables, hfd_tables, 11)
water_df <- na.omit(water_df)

df_water_barplot_cd <- transform_for_barplot(water_df,
                                             "water_consumption",
                                             "CD",
                                             1)

df_water_barplot_cd_mean <- df_water_barplot_cd %>%
  group_by(phase, genotype) %>%
  summarise(mean = mean(water_consumption),
            sd = sd(water_consumption))



bodycomp_df <- merge_cd_hfd(cd_tables, hfd_tables, 13)
bodycomp_df <- na.omit(bodycomp_df)

df_bodycomp_barplot_cd <- transform_for_barplot(bodycomp_df,
                                                "mass", "CD", 2)

df_bodycomp_barplot_cd_mean <- df_bodycomp_barplot_cd %>%
  group_by(phase, genotype) %>%
  summarise(mean = mean(mass),
            sd = sd(mass))



bodycomp_norm_df <- merge_cd_hfd(cd_tables, hfd_tables, 14)
bodycomp_norm_df <- na.omit(bodycomp_norm_df)

df_bodycomp_norm_barplot_cd <- transform_for_barplot(bodycomp_norm_df,
                                                     "mass", "CD", 2)

df_bodycomp_norm_barplot_cd_mean <- df_bodycomp_norm_barplot_cd %>%
  group_by(phase, genotype) %>%
  summarise(mean = mean(mass),
            sd = sd(mass))

activity_df <- merge_cd_hfd(cd_tables, hfd_tables, 8)
activity_df <- na.omit(activity_df)

df_activity_barplot_cd <- transform_for_barplot(activity_df,
                                                "breaks", "CD", 1)

df_activity_barplot_cd_mean <- df_activity_barplot_cd %>%
  group_by(phase, genotype) %>%
  summarise(mean = mean(breaks),
            sd = sd(breaks))

zrearing_df <- merge_cd_hfd(cd_tables, hfd_tables, 9)
zrearing_df <- na.omit(zrearing_df)

df_zrearing_barplot_cd <- transform_for_barplot(zrearing_df,
                                                "breaks", "CD", 1)

df_zrearing_barplot_cd_mean <- df_zrearing_barplot_cd %>%
  group_by(phase, genotype) %>%
  summarise(mean = mean(breaks),
            sd = sd(breaks))

ebalance_df <- merge_cd_hfd(cd_tables, hfd_tables, 15)
ebalance_df <- na.omit(ebalance_df)

df_ebalance_barplot_cd <- transform_for_barplot(ebalance_df,
                                                "kcal", "CD", 1)

df_ebalance_barplot_cd_mean <- df_ebalance_barplot_cd %>%
  group_by(phase, genotype) %>%
  summarise(mean = mean(kcal),
            sd = sd(kcal))


EE_phases_df <- merge_cd_hfd(cd_tables, hfd_tables, 7)
EE_phases_df <- na.omit(EE_phases_df)

df_EE_phases_barplot_cd <- transform_for_barplot(EE_phases_df,
                                                "kcal", "CD", 1)

df_EE_phases_barplot_cd_mean <- df_EE_phases_barplot_cd %>%
  group_by(phase, genotype) %>%
  summarise(mean = mean(kcal),
            sd = sd(kcal))


food_phase_df <- merge_cd_hfd(cd_tables, hfd_tables, 10)
food_phase_df <- na.omit(food_phase_df)

df_food_phase_barplot_cd <- transform_for_barplot(food_phase_df,
                                                 "kcal", "CD", 1)

df_food_phase_barplot_cd_mean <- df_food_phase_barplot_cd %>%
  group_by(phase, genotype) %>%
  summarise(mean = mean(kcal),
            sd = sd(kcal))

#### RER ------
rer_phase_df <- merge_cd_hfd(cd_tables, hfd_tables, 6)
rer_phase_df <- na.omit(rer_phase_df)

df_rer_phase_barplot_cd <- transform_for_barplot(rer_phase_df,
                                                 "rq", "CD", 1)

df_rer_phase_barplot_cd_mean <- df_rer_phase_barplot_cd %>%
  group_by(phase, genotype) %>%
  summarise(mean = mean(rq),
            sd = sd(rq))

### HFD -----

df_water_barplot_hfd <- transform_for_barplot(water_df,
                                             "water_consumption", "HFD", 1)

df_water_barplot_hfd_mean <- df_water_barplot_hfd %>%
  group_by(phase, genotype) %>%
  summarise(mean = mean(water_consumption),
            sd = sd(water_consumption))



df_bodycomp_barplot_hfd <- transform_for_barplot(bodycomp_df,
                                                "mass", "HFD", 2)

df_bodycomp_barplot_hfd_mean <- df_bodycomp_barplot_hfd %>%
  group_by(phase, genotype) %>%
  summarise(mean = mean(mass),
            sd = sd(mass))



df_bodycomp_norm_barplot_hfd <- transform_for_barplot(bodycomp_norm_df,
                                                     "mass", "HFD", 2)

df_bodycomp_norm_barplot_hfd_mean <- df_bodycomp_norm_barplot_hfd %>%
  group_by(phase, genotype) %>%
  summarise(mean = mean(mass),
            sd = sd(mass))


df_activity_barplot_hfd <- transform_for_barplot(activity_df,
                                                "breaks", "HFD", 1)

df_activity_barplot_hfd_mean <- df_activity_barplot_hfd %>%
  group_by(phase, genotype) %>%
  summarise(mean = mean(breaks),
            sd = sd(breaks))

df_zrearing_barplot_hfd <- transform_for_barplot(zrearing_df,
                                                "breaks", "HFD", 1)

df_zrearing_barplot_hfd_mean <- df_zrearing_barplot_hfd %>%
  group_by(phase, genotype) %>%
  summarise(mean = mean(breaks),
            sd = sd(breaks))

df_ebalance_barplot_hfd <- transform_for_barplot(ebalance_df,
                                                "kcal", "HFD", 1)

df_ebalance_barplot_hfd_mean <- df_ebalance_barplot_hfd %>%
  group_by(phase, genotype) %>%
  summarise(mean = mean(kcal),
            sd = sd(kcal))


df_EE_phases_barplot_hfd <- transform_for_barplot(EE_phases_df,
                                                 "kcal", "HFD", 1)

df_EE_phases_barplot_hfd_mean <- df_EE_phases_barplot_hfd %>%
  group_by(phase, genotype) %>%
  summarise(mean = mean(kcal),
            sd = sd(kcal))


df_food_phase_barplot_hfd <- transform_for_barplot(food_phase_df,
                                                  "kcal", "HFD", 1)

df_food_phase_barplot_hfd_mean <- df_food_phase_barplot_hfd %>%
  group_by(phase, genotype) %>%
  summarise(mean = mean(kcal),
            sd = sd(kcal))

#### RER ------
df_rer_phase_barplot_hfd <- transform_for_barplot(rer_phase_df,
                                                 "rq", "HFD", 1)

df_rer_phase_barplot_hfd_mean <- df_rer_phase_barplot_hfd %>%
  group_by(phase, genotype) %>%
  summarise(mean = mean(rq),
            sd = sd(rq))

# Plots -------------------------------------------------------------------
setwd("/Users/artemii/Yandex.Disk.localized/Inserm/Exps/metabolic-cages-feb-2021/metabolic-cages-2021/plots")

## Plot BW vs EE 24h ------------

EE_df_cd_wt <- EE_df[which(EE_df$genotype=='WT' & EE_df$diet=='CD'), ]
EE_df_cd_ko <- EE_df[which(EE_df$genotype=='KO' & EE_df$diet=='CD'), ]
EE_df_hfd_wt <- EE_df[which(EE_df$genotype=='WT' & EE_df$diet=='HFD'), ]
EE_df_hfd_ko <- EE_df[which(EE_df$genotype=='KO' & EE_df$diet=='HFD'), ]
EE_df_cd <- EE_df[which(EE_df$diet=='CD'), ]
EE_df_hfd <- EE_df[which(EE_df$diet=='HFD'), ]

ee_size = 5
ee_stroke = 1.5

### HFD -----
ggplot() +
  #geom_point(data = EE_df_cd_wt,
  #           mapping = aes(x = BW, y = EE_24h),
  #           colour = 'black') +
  #geom_point(data = EE_df_cd_ko,
  #           mapping = aes(x = BW, y = EE_24h),
  #           colour = 'red') +
  geom_smooth(data = EE_df_cd_wt,
              mapping = aes(x = BW, y = EE_24h),
              method = "lm",
              se = TRUE,
              colour = '#B9B9B9') +
  geom_smooth(data = EE_df_cd_ko,
              mapping = aes(x = BW, y = EE_24h),
              method = "lm",
              se = TRUE,
              colour = '#009193') +
  geom_point(data = EE_df_cd_wt,
             mapping = aes(x = BW, y = EE_24h),
             colour = '#B9B9B9', size = ee_size) +
  geom_point(data = EE_df_cd_wt,
             mapping = aes(x = BW, y = EE_24h),
             shape = 1,
             colour = 'black', size = ee_size, stroke = ee_stroke) +
  geom_point(data = EE_df_cd_ko,
             mapping = aes(x = BW, y = EE_24h),
             colour = '#009193', size = ee_size) +
  geom_point(data = EE_df_cd_ko,
             mapping = aes(x = BW, y = EE_24h),
             shape = 1,
             colour = 'black', size = ee_size, stroke = ee_stroke) +
  theme_classic(base_size = 20) +
  scale_y_continuous(limits = c(8,17.5))
#  scale_x_continuous(limits = c(20,50))

ggsave("ee_vs_bw_cd.svg", width = 10, height = 10, units = "cm")
dev.off()

### CD -----

ggplot() +
  #geom_point(data = EE_df_cd_wt,
  #           mapping = aes(x = BW, y = EE_24h),
  #           colour = 'black') +
  #geom_point(data = EE_df_cd_ko,
  #           mapping = aes(x = BW, y = EE_24h),
  #           colour = 'red') +
  geom_smooth(data = EE_df_hfd_wt,
              mapping = aes(x = BW, y = EE_24h),
              method = "lm",
              se = TRUE,
              colour = '#B9B9B9') +
  geom_smooth(data = EE_df_hfd_ko,
              mapping = aes(x = BW, y = EE_24h),
              method = "lm",
              se = TRUE,
              colour = '#009193') +
  geom_point(data = EE_df_hfd_wt,
             mapping = aes(x = BW, y = EE_24h),
             colour = '#B9B9B9', size = ee_size) +
  geom_point(data = EE_df_hfd_wt,
             mapping = aes(x = BW, y = EE_24h),
             shape = 1,
             colour = 'black', size = ee_size, stroke = ee_stroke) +
  geom_point(data = EE_df_hfd_ko,
             mapping = aes(x = BW, y = EE_24h),
             colour = '#009193', size = ee_size) +
  geom_point(data = EE_df_hfd_ko,
             mapping = aes(x = BW, y = EE_24h),
             shape = 1,
             colour = 'black', size = ee_size, stroke = ee_stroke) +
  theme_classic(base_size = 20) +
  scale_y_continuous(limits = c(8,17.5))
  
ggsave("ee_vs_bw_hfd.svg", width = 10, height = 10, units = "cm")
dev.off()
  #geom_point(data = EE_df[(EE_df$genotype == 'KO') && (EE_df$diet == 'CD'),],
  #           mapping = aes(x = BW, y = EE_24h),
  #           colour = 'red', size = 3)



## Plotting time series -------

plot_time_series(RER_df, diet = "cd", use_scaling = T,
                 min_y_val = 0.65, max_y_val = 1)
ggsave("rer_cd_series.svg", width = 14, height = 10, units = "cm")
dev.off()

plot_time_series(RER_df, diet = "hfd", use_scaling = T,
                 min_y_val = 0.65, max_y_val = 1)
ggsave("rer_hfd_series.svg", width = 14, height = 10, units = "cm")
dev.off()

plot_time_series(EE_df, diet = "cd", use_scaling = T,
                 min_y_val = 0.3, max_y_val = 0.8)
ggsave("ee_cd_series.svg", width = 14, height = 10, units = "cm")
dev.off()


plot_time_series(EE_df, diet = "hfd", use_scaling = T,
                 min_y_val = 0.3, max_y_val = 0.8)
ggsave("ee_hfd_series.svg", width = 14, height = 10, units = "cm")
dev.off()


plot_time_series(food, diet = "cd", use_scaling = T,
                 min_y_val = 0, max_y_val = 18)
ggsave("food_cons_cd_series.svg", width = 14, height = 10, units = "cm")
dev.off()

plot_time_series(food, diet = "hfd", use_scaling = T,
                 min_y_val = 0, max_y_val = 18)
ggsave("food_cons_hfd_series.svg", width = 14, height = 10, units = "cm")
dev.off()

## Plotting bar plots -------
### Weight ----------

ggplot(data = BW_mean_barplot, aes(x = diet_mean, y = mean_bw,
                                   group = genotype_vector_mean)) +
  geom_col(aes(x = diet_mean),
           position = position_dodge(), colour = "black",
           fill = "white",
           size = 0.8) +
  geom_jitter(data = BW_df_barplot, aes(x = diet, y = BW,
                                        fill= genotype_vector),
             na.rm = T,
             position = position_jitterdodge(dodge.width = 1,
                                             jitter.width = 0.3),
             shape = 21, size = 3,
             colour = "black", stroke = 1.5,
             inherit.aes = F) +
  scale_fill_manual(values = c(WT = "#B9B9B9", KO = "#009193")) +
  geom_errorbar(aes(ymin = mean_bw - se,
                    ymax = mean_bw + se),
                position = position_dodge(0.9),
                colour = "grey",
                alpha = 1,
                width = 0.3,
                size = 0.5) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 50)) +
  scale_x_discrete(expand = c(0,0))

ggsave("body_weight.svg", width = 6.23, height = 10, units = "cm")
dev.off()

### Phase plots (Water consumption and alike) ------


#### CD plots ------
plot_phase_barplot(df_water_barplot_cd,
                   df_water_barplot_cd_mean,
                   var_num = 1,
                   scale_factor = 0,
                   min_val = 0,
                   max_val = 0.31)
ggsave("water_cons_cd_phases.svg", width = 8, height = 10, units = "cm")
dev.off()


plot_phase_barplot(df_food_phase_barplot_cd,
                   df_food_phase_barplot_cd_mean,
                   var_num = 4,
                   scale_factor = 0.65)
ggsave("food_cons_cd_phases.svg", width = 8, height = 10, units = "cm")
dev.off()

plot_phase_barplot(df_activity_barplot_cd,
                   df_activity_barplot_cd_mean,
                   var_num = 3,
                   scale_factor = 0,
                   min_val = 0,
                   max_val = 1600)
ggsave("activity_cd_phases.svg", width = 8, height = 10, units = "cm")
dev.off()


plot_phase_barplot(df_zrearing_barplot_cd,
                   df_zrearing_barplot_cd_mean,
                   var_num = 3,
                   scale_factor = 0.15)
ggsave("zrearing_cd_phases.svg", width = 8, height = 10, units = "cm")
dev.off()

##### Body composition -----
plot_phase_barplot(df_bodycomp_barplot_cd,
                   df_bodycomp_barplot_cd_mean,
                   var_num = 2,
                   scale_factor = 0,
                   min_val = 0,
                   max_val = 70)
ggsave("bodycomp_cd_phases.svg", width = 8, height = 10, units = "cm")
dev.off()


##### Normalized Body composition  -----
plot_phase_barplot(df_bodycomp_norm_barplot_cd,
                   df_bodycomp_norm_barplot_cd_mean,
                   var_num = 2,
                   scale_factor = 0,
                   min_val = 0,
                   max_val = 80)
ggsave("bodycomp_norm_cd_phases.svg", width = 8, height = 10, units = "cm")
dev.off()


plot_phase_barplot(df_EE_phases_barplot_cd,
                   df_EE_phases_barplot_cd_mean,
                   var_num = 4,
                   scale_factor = 1.5)
ggsave("EE_cd_phases.svg", width = 8, height = 10, units = "cm")
dev.off()


plot_phase_barplot(df_ebalance_barplot_cd,
                   df_ebalance_barplot_cd_mean,
                   var_num = 4,
                   scale_factor = 0,
                   min_val = -0.5,
                   max_val = 0.6)
ggsave("ebalance_cd_phases.svg", width = 8, height = 10, units = "cm")
dev.off()

##### RER  -----
# I don't understand why it doesn't work. I'll just draw the bars myself
plot_phase_barplot(df_rer_phase_barplot_cd,
                   df_rer_phase_barplot_cd_mean,
                   var_num = 5,
                   scale_factor = 0,
                   min_val = 0.6,
                   max_val = 1)
ggsave("rer_cd_phases.svg", width = 8, height = 10, units = "cm")
dev.off()

#### HFD plots ------

plot_phase_barplot(df_water_barplot_hfd,
                   df_water_barplot_hfd_mean,
                   var_num = 1,
                   scale_factor = 0,
                   min_val = 0,
                   max_val = 0.31)
ggsave("water_cons_hfd_phases.svg", width = 8, height = 10, units = "cm")
dev.off()


plot_phase_barplot(df_food_phase_barplot_hfd,
                   df_food_phase_barplot_hfd_mean,
                   var_num = 4,
                   scale_factor = 0.65)
ggsave("food_cons_hfd_phases.svg", width = 8, height = 10, units = "cm")
dev.off()


plot_phase_barplot(df_activity_barplot_hfd,
                   df_activity_barplot_hfd_mean,
                   var_num = 3,
                   scale_factor = 0,
                   min_val = 0,
                   max_val = 1600)
ggsave("activity_hfd_phases.svg", width = 8, height = 10, units = "cm")
dev.off()



plot_phase_barplot(df_zrearing_barplot_hfd,
                   df_zrearing_barplot_hfd_mean,
                   var_num = 3,
                   scale_factor = 0.15)
ggsave("zrearing_hfd_phases.svg", width = 8, height = 10, units = "cm")
dev.off()


##### Body composition -----
plot_phase_barplot(df_bodycomp_barplot_hfd,
                   df_bodycomp_barplot_hfd_mean,
                   var_num = 2,
                   scale_factor = 0,
                   min_val = 0,
                   max_val = 70)
ggsave("bodycomp_hfd_phases.svg", width = 8, height = 10, units = "cm")
dev.off()


##### Normalized Body composition  -----
plot_phase_barplot(df_bodycomp_norm_barplot_hfd,
                   df_bodycomp_norm_barplot_hfd_mean,
                   var_num = 2,
                   scale_factor = 0,
                   min_val = 0,
                   max_val = 80)
ggsave("bodycomp_norm_hfd_phases.svg", width = 8, height = 10, units = "cm")
dev.off()



plot_phase_barplot(df_EE_phases_barplot_hfd,
                   df_EE_phases_barplot_hfd_mean,
                   var_num = 4,
                   scale_factor = 1.5)
ggsave("EE_hfd_phases.svg", width = 8, height = 10, units = "cm")
dev.off()


plot_phase_barplot(df_ebalance_barplot_hfd,
                   df_ebalance_barplot_hfd_mean,
                   var_num = 4,
                   scale_factor = 0,
                   min_val = -0.5,
                   max_val = 0.6)
ggsave("ebalance_hfd_phases.svg", width = 8, height = 10, units = "cm")
dev.off()

##### RER  -----

plot_phase_barplot(df_rer_phase_barplot_hfd,
                   df_rer_phase_barplot_hfd_mean,
                   var_num = 5,
                   scale_factor = 0,
                   min_val = 0.6,
                   max_val = 1)
ggsave("rer_hfd_phases.svg", width = 8, height = 10, units = "cm")
dev.off()

# Statistics --------------------------------------------------------------

# ANCOVA tests for body weight normalized energy expenditure
library(multcomp)

ancova_hfd <- aov(EE_24h ~ genotype*BW, data = EE_df_hfd)
summary(ancova_hfd)

tukey.test <- TukeyHSD(ancova_hfd, which = "genotype")
tukey.test

ancova_cd <- aov(EE_24h ~ genotype*BW, data = EE_df_cd)
summary(ancova_cd)

tukey.test <- TukeyHSD(ancova_cd, which = "genotype")
tukey.test

# ANOVA with repeated measures for RER


