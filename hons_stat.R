#HOW DOES THE POTATO HYPHOSPHERE PHOSPHORUS AVAILABILITY RESPOND TO DROUGHT STRESS
#Author: Hazel Surtees
#Contact: h.e.surtees@dundee.ac.uk
#Last Updated: 17.01.24

#Data is available from Prof. Tim George, tim.george@hutton.ac.uk

##Load required packages
library(ggplot2)
library(readxl)
library(dplyr)
library(ggpattern)
library(ggpubr)

##load and attach data

DataALL3 <- read_excel("DataALL3.xlsx", sheet = "Sheet2")
View(DataALL3)
attach(DataALL3)

r_yield <- read_excel("r_yield.xlsx")
View(r_yield)
attach(r_yield)

amfper <- read_excel("DataALL3.xlsx", sheet = "Sheet1")
View(amfper)
attach(amfper)

platesALL <- read_excel("plate results/hazel_platesALL.xlsx", 
                        sheet = "Po")
View(platesALL)
attach(platesALL)

##combine data.
### note that most of this data relates directly to soil samples, n = 39,
### however the phoD analysis was conducted on DNA extracted in triplicate from the soils samples
### hence analysis of references a slightly differently organised datafram
### doubtless, there's a more elegant solution. but this works.

soil.df = subset(DataALL3, select = -c(total, Po, olsen, phosphatase))
View(soil.df)
phos_join <- inner_join(platesALL, soil.df, by = c("field_pos", "amf_innoc"), multiple = "all")
head(phos_join)

##check for normal dist
hist(DataALL3$soil_mois)
hist(DataALL3$pH)
hist(phos_join$Po)
hist(phos_join$Pi)
hist(phos_join$Tot_P)
hist(phos_join$phoD)

## filter nonsense values, unaligned data (NAs), and transform
##tranforming phoD by log10 and converting phosphatase activity to katal units.
amfper <- amfper %>%
  filter(water != "NA") %>%
  mutate(treatment = case_when(             
    grepl("Drought", water) ~ "Drought",
    grepl("Irr", water) ~ "Irrigated",
  )
  )

phos_join.g <- phos_join %>%   # use this for tests on soil samples
  mutate(phtase_katal = phosphatase/(139.11*3.6)) %>%
  filter(sample == "Soil") %>%
  filter(primer =="16S") %>%
  filter(Po > 0) %>%
  filter(Pi > 0) %>%
  filter(phosphatase > 0) %>%
  mutate(treatment = case_when(               # creates the genus column and specifies conditions
    grepl("drought", water) & grepl("P", amf_innoc) ~ "AMF & Drought",
    grepl("drought", water) & grepl("Con", amf_innoc) ~ "Control & Drought",
    grepl("irr", water) & grepl("P", amf_innoc) ~ "AMF & Irrigated",
    grepl("irr", water) & grepl("Con", amf_innoc) ~ "Control & Irrigated",
  )
  )

phod.g <- phos_join %>%  #use this for phoD tests which can be linked to 96well
  mutate(phtase_katal = phosphatase/(139.11*3.6)) %>%
  filter(Tot_P > 0) %>%
  filter(Po > 0) %>%
  filter(Pi > 0) %>%
  filter(phoD > 0) %>%
  mutate(log_phoD = log(phoD)) %>%
  filter(phosphatase > 0) %>%
  mutate(treatment = case_when(               # creates the genus column and specifies conditions
    grepl("drought", water) & grepl("P", amf_innoc) ~ "AMF & Drought",
    grepl("drought", water) & grepl("Con", amf_innoc) ~ "Control & Drought",
    grepl("irr", water) & grepl("P", amf_innoc) ~ "AMF & Irrigated",
    grepl("irr", water) & grepl("Con", amf_innoc) ~ "Control & Irrigated",
  )
  )

##check how many values filtered out
preprune_soil <- length(phos_join$sample_ID)
postprune_soil <- length(phod.g$sample_ID)
postprune_phod <- length(phos_join.g$sample_ID)


##statistical analysis - replace x & y as necessary (I'll learn loops one day)
x <- phod.g$log_phoD
y <- phos_join.g$phtase_katal

anova_totP<- aov(y ~ treatment, data = phos_join.g)
anova_Po<- aov(y ~ treatment, data = phos_join.g)
anova_Pi<- aov(y ~ treatment, data = phos_join.g)
anova_phtase<- aov(y ~ treatment, data = phos_join.g)
anova_phod <- aov(x ~ treatment, data = phod.g)
anova_pH <- aov(y ~ treatment, data = phos_join.g)
anova_mois <- aov(y ~ treatment, data = phos_join.g)
t_peramf <- t.test(per_amf ~ treatment, data = amfper)

t_yield <- t.test(yield ~ treatment, data = r_yield)
anova_yield <- aov(yield ~ treatment + cultivar + treatment:cultivar, data = r_yield)


summary(anova_totP)
summary(anova_Po)
summary(anova_Pi)
summary(anova_phtase)
summary(anova_pH)
summary(anova_mois)
summary(anova_phod)
t_peramf
summary(anova_yield)


#tukey hsd

TukeyHSD(anova_totP)
TukeyHSD(anova_Pi)
TukeyHSD(anova_phtase)
TukeyHSD(anova_pH)
TukeyHSD(anova_mois)

#effect size (treatment sum of squares/(treatment sum of sq + res sum of sq)
es_totP <- 484.8/(484.8+1403.5)
es_Pi <- 319.6/(319.6+1059.0)
es_phtase <- 0.004491/(0.004491+0.014140)
es_pH <- 0.03171/(0.03171+0.14573)
es_mois <- 825.1/(825.1+2746.2)

##Plotting

(yield_box.ggstat <-
    ggplot(r_yield, aes(x = treatment, y = yield, fill = treatment)) +
    geom_boxplot(aes(fill = treatment),linewidth = 0.7)+
    scale_fill_manual(name = "treatment", values = c("#F4A460","#63B8FF")) +
    stat_compare_means(method = "t.test", label.y =9) +
    stat_compare_means(label = "p.signif", method = "t.test", ref.group =".all.", hide.ns = TRUE) +
    geom_hline(yintercept = mean(r_yield$yield), linetype = 2)+
    theme_bw() +
    scale_y_continuous(limits = c(0, 10)) +
    ylab("Total Yield (Kg)") +
    theme(axis.title.y = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 10, face = "bold", colour = "black"),
          axis.title.x=element_blank(), 
          plot.margin = unit(c(1,1,1,1), units = , "cm")))

#Test Drought Effect on AMF-root colonisation

(per_amf_box.ggstat <-
    ggplot(amfper, aes(x = treatment, y = per_amf, fill = treatment)) +
    geom_boxplot(aes(fill = treatment),linewidth = 1)+
    scale_fill_manual(name = c("Drought", "Irrigated"), values = c("#F4A460","#63B8FF")) +
    stat_compare_means(method = "t.test", label.y =27) +
    stat_compare_means(label = "p.signif", method = "t.test", ref.group =".all.", hide.ns = TRUE) +
    geom_hline(yintercept = mean(amfper$per_amf), linetype = 2)+
    labs(x = (name = c("Drought", "Irrigated"))) +
    theme_bw() +
    scale_y_continuous(limits = c(0, 30)) +
    ylab("% AMF Colonisation") +
    theme(axis.title.y = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 12, face = "bold", colour = "black"),
          axis.text.x = element_text(),
          axis.title.x = element_blank(),
          plot.margin = unit(c(1,1,1,1), units = , "cm")))

##Test AMF/drought effect on pH
(pH.ggstat <-
    ggplot(phos_join.g, aes(x = treatment, y = pH, pattern = amf_innoc, fill = water)) +
    geom_boxplot(aes(fill = water),linewidth = 1)+
    scale_fill_manual(name = "water", values = c("#F4A460","#63B8FF")) +
    stat_compare_means(method = "anova", label.y =5.45) +
    stat_compare_means(label = "p.signif", method = "t.test", ref.group =".all.", hide.ns = TRUE) +
    geom_hline(yintercept = mean(phos_join.g$pH), linetype = 2)+
    geom_boxplot_pattern(position = position_dodge(preserve = "single"), color = "black", pattern_fill = "white", pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
    guides(pattern = guide_legend(override.aes = list(fill = "white")), fill = guide_legend(override.aes = list(pattern = "none"))) +
    theme_bw() +
    scale_y_continuous(limits = c(5.1, 5.5)) +
    ylab("Soil pH") +
    theme(axis.title.y = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 7, face = "bold", colour = "black"),
          axis.title.x=element_blank(), 
          plot.margin = unit(c(1,1,1,1), units = , "cm")))

##Test AMF/drought effect on soil moisture
(moisture.ggstat <-
    ggplot(phos_join.g, aes(x = treatment, y = soil_mois, pattern = amf_innoc, fill = water)) +
    geom_boxplot(aes(fill = water),linewidth = 1)+
    scale_fill_manual(name = "water", values = c("#F4A460","#63B8FF")) +
    stat_compare_means(method = "anova", label.y =56) +
    stat_compare_means(label = "p.signif", method = "t.test", ref.group =".all.", hide.ns = TRUE) +
    geom_hline(yintercept = mean(phos_join.g$soil_mois), linetype = 2)+
    geom_boxplot_pattern(position = position_dodge(preserve = "single"), color = "black", pattern_fill = "white", pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
    guides(pattern = guide_legend(override.aes = list(fill = "white")), fill = guide_legend(override.aes = list(pattern = "none"))) +
    theme_bw() +
    scale_y_continuous(limits = c(0, 60)) +
    ylab("Soil Moisture (Weight Difference (g))") +
    theme(axis.title.y = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 7, face = "bold", colour = "black"),
          axis.title.x=element_blank(), 
          plot.margin = unit(c(1,1,1,1), units = , "cm")))

##Test AMF/drought effect on soil phosphorus and phod
(tot_p_box.ggstat <-
    ggplot(phos_join.g, aes(x = treatment, y = Tot_P, pattern = amf_innoc, fill = water)) +
    geom_boxplot(aes(fill = water),linewidth = 1)+
    scale_fill_manual(name = "water", values = c("#F4A460","#63B8FF")) +
    stat_compare_means(method = "anova", label.y =40) +
    stat_compare_means(label = "p.signif", method = "t.test", ref.group =".all.", hide.ns = TRUE) +
    geom_hline(yintercept = mean(phos_join.g$Tot_P), linetype = 2)+
    geom_boxplot_pattern(position = position_dodge(preserve = "single"), color = "black", pattern_fill = "white", pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
    guides(pattern = guide_legend(override.aes = list(fill = "white")), fill = guide_legend(override.aes = list(pattern = "none"))) +
    theme_bw() +
    scale_y_continuous(limits = c(0, 45)) +
    ylab("Total Phosphorus (μg/g soil)") +
    theme(axis.title.y = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 7, face = "bold", colour = "black"),
          axis.title.x=element_blank(), 
          plot.margin = unit(c(1,1,1,1), units = , "cm")))


(po_box.ggstat <-
    ggplot(phos_join.g, aes(x = treatment, y = Po, pattern = amf_innoc, fill = water)) +
    geom_boxplot(aes(fill = water),linewidth = 1)+
    scale_fill_manual(name = "water", values = c("#F4A460","#63B8FF")) +
    stat_compare_means(method = "anova", label.y =9) +
    stat_compare_means(label = "p.signif", method = "t.test", ref.group =".all.", hide.ns = TRUE) +
    geom_boxplot_pattern(position = position_dodge(preserve = "single"), color = "black", pattern_fill = "white", pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
    guides(pattern = guide_legend(override.aes = list(fill = "white")), fill = guide_legend(override.aes = list(pattern = "none"))) +
    geom_hline(yintercept = mean(phos_join.g$Po), linetype = 2)+
    theme_bw() +
    scale_y_continuous(limits = c(0, 10)) +
    ylab("Organic Phosphorus (μg/g soil)") +
    theme(axis.title.y = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 7, face = "bold", colour = "black"),
          axis.title.x=element_blank(), 
          plot.margin = unit(c(1,1,1,1), units = , "cm")))

(pi_box.ggstat <-
    ggplot(phos_join.g, aes(x = treatment, y = Pi, pattern = amf_innoc, fill = water)) +
    geom_boxplot(aes(fill = water),linewidth = 1)+
    scale_fill_manual(name = "water", values = c("#F4A460","#63B8FF")) +
    stat_compare_means(method = "anova", label.y =30) +
    stat_compare_means(label = "p.signif", method = "t.test", ref.group =".all.", hide.ns = TRUE) +
    geom_boxplot_pattern(position = position_dodge(preserve = "single"), color = "black", pattern_fill = "white", pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
    guides(pattern = guide_legend(override.aes = list(fill = "white")), fill = guide_legend(override.aes = list(pattern = "none"))) +
    geom_hline(yintercept = mean(phos_join.g$Pi), linetype = 2)+
    theme_bw() +
    #    scale_y_continuous(limits = c(11.9, 12.6), breaks = seq(11.9, 12.6, 0.1)) +
    ylab("Inorganic Phosphorus (μg/g soil)") +
    theme(axis.title.y = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 7, face = "bold", colour = "black"),
          axis.title.x=element_blank(), 
          plot.margin = unit(c(1,1,1,1), units = , "cm")))

(phtase_box.gg <-
    ggplot(phos_join.g, aes(x = treatment, y = phtase_katal, pattern = amf_innoc, fill = water)) +
    geom_boxplot(aes(fill = water),linewidth = 1)+
    scale_fill_manual(name = "water", values = c("#F4A460","#63B8FF")) +
    stat_compare_means(method = "anova", label.y =0.1) +
    stat_compare_means(label = "p.signif", method = "t.test", ref.group =".all.", hide.ns = TRUE) +
    geom_boxplot_pattern(position = position_dodge(preserve = "single"), color = "black", pattern_fill = "white", pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
    guides(pattern = guide_legend(override.aes = list(fill = "white")), fill = guide_legend(override.aes = list(pattern = "none"))) +
    geom_hline(yintercept = mean(phos_join.g$phtase_katal), linetype = 2)+
    theme_bw() +
    #    scale_y_continuous(limits = c(11.9, 12.6), breaks = seq(11.9, 12.6, 0.1)) +
    ylab(bquote(bold('Phosphatase Activity'~(nKatg^-1)))) +
    theme(axis.title.y = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 7, face = "bold", colour = "black"),
          axis.title.x=element_blank(), 
          plot.margin = unit(c(1,1,1,1), units = , "cm")))

(phod_boxplot <-
    ggplot(phod.g, aes(x = treatment, y = log_phoD, pattern = amf_innoc, fill = water)) +
    geom_boxplot(aes(fill = water), linewidth = 1)+
    scale_fill_manual(name = "water", values = c("#F4A460","#63B8FF")) +
    stat_compare_means(method = "anova", label.y =12.6) +
    geom_boxplot_pattern(position = position_dodge(preserve = "single"), color = "black", pattern_fill = "white", pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
    guides(pattern = guide_legend(override.aes = list(fill = "white")), fill = guide_legend(override.aes = list(pattern = "none"))) +
    geom_hline(yintercept = mean(phod.g$log_phoD), linetype = 2)+
    theme_bw() +
    scale_y_continuous(limits = c(11.9, 12.6), breaks = seq(11.9, 12.6, 0.1)) +
    ylab("Log Transformed PhoD Concentration") +
    theme(axis.title.y = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 7, face = "bold", colour = "black"),
          axis.title.x=element_blank(), 
          plot.margin = unit(c(1,1,1,1), units = , "cm")))   
