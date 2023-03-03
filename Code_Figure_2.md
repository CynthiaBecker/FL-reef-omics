---
title: "Code to Reproduce Figure 2"
subtitle: "Microorganisms and dissolved metabolites distinguish Florida’s Coral Reef habitats"
author: "Cynthia Becker"
date: "2023-03-03"
output:
  html_document:
    keep_md: true
    toc: true
    toc_depth: 3
    number_sections: true
    theme: united
---



# Setup
### Install necessary packages


```r
# For data wrangling
library(dplyr); packageVersion("dplyr")
library(tidyr); packageVersion("tidyr")

# For microbiome/metagenome analysis

# For visualization
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())
```

### Read in prepped data

This table is required for the subsequent analyses. Load prior to generating figures.


```r
data2 <- read.delim("data/FLK2019_environmental_data_for_R.txt", header = TRUE, sep = "\t") 

# Fix the data frame a bit so zone is a factor and no3 and silicate are numeric
data2 <- as.data.frame(data2)
data2$zonenumber <- as.factor(data2$zonenumber)
data2$no3 <- as.numeric(data2$no3)
data2$silicate <- as.numeric(data2$silicate)
```


# Figure 2

**Zones in Florida’s Coral Reef (FCR) have significantly different disease prevalence, microbial  abundances and organic carbon and silicate concentrations.**

**Figure 2.** a) Map of all 85 reefs sampled for environmental variables are colored by zone and placed above FCR (tan). Environmental variables that changed significantly by zone include b) stony coral tissue loss disease prevalence (percent of live coral), (c) ratio of heterotrophic to photosynthetic microbial cells, (d) total organic carbon concentration, (e) heterotrophic microbial abundances (unpigmented bacteria and archaea), (f) Prochlorococcus abundances, (g) Synechococcus abundances, (h) picoeukaryote abundances, and (i) silicate concentration was significantly different between FCR zones as a result of a Kruskal-Wallis test and less than the Bonferroni-corrected p-value of 0.00192. Asterisks denote the zone was significantly different than all other zones tested in the pairwise Wilcoxon rank sum test (Benjamini-Hochberg adjusted p<0.05). 


```r
# I just produced boxplots of the Environmental parameters across each environemntal zone. I will then conduct a Kruskal-Wallis test to test for significant differences between the values in each Zone or Region.
# About the Kruskal-Wallis test: Kruskal-Wallis test by rank is a non-parametric alternative to one-way ANOVA test, which extends the two-samples Wilcoxon test in the situation where there are more than two groups (I have several groups since I am testing for zone differences). It’s recommended when the assumptions of one-way ANOVA test are not met. These assumptions include equal variance and normal distribution at each factor level. These are not too bad if not true, but if sample sizes are not equal that also makes the ANOVA test more challenging. 

data3 <- data2 %>%
  select(-temp, -dissolvedO2_percent, -dissolvedO2, -salinity, -Metagenomes, -Microbiomes, -CoralTissue, -MicrobiomeANDTissue, -BenthicSurvey) %>%
  mutate(CellRatio = hbact / (pro + syn)) %>% #get ratio of heterotrophs to autotrophic bacteria
  mutate(TON = tn - no2no3 - nh4) %>% #get total organic nitrogen
  mutate(All_algae = Calcified_Macroalgae + Fleshy_Macroalgae + Turf_Algae) %>%
  mutate(All_macroalgae = Calcified_Macroalgae + Fleshy_Macroalgae) %>%
  mutate(MA_Turf = Fleshy_Macroalgae + Turf_Algae) %>%
  mutate(Coral_to_Algae = Hard_coral/MA_Turf) %>%
  select(-no2no3)

#Goal: Loop tyhrough each environmental variable. Do a kruskal-wallis test on each of the variables. Save the p-value. Also generate a bonferroni-corrected p-value because you are doing multiple comparisons and the likelihood of getting a significant results due to chance is high. If the bonferroni-corrected p-value is less than 0.05, conduct a post test that is a pairwise wilcoxon test to get p-values between each combination of zones. 

variables <- colnames(data3[,18:50])
kruskal.result <- data.frame()
wilcox.result <- list()
Bonf.alpha <- 0.05/length(variables)

for (i in 18:50) {
  result <- kruskal.test(data3[,i] ~ data3[,6]) #do kruskal wallis test on each variable of interest
  kruskal.result <- rbind(kruskal.result, result) #save kruskal wallis test results
  wilcox <- pairwise.wilcox.test(data3[ ,i], data3[ ,6], p.adjust.method = "BH") #use pairwise test to see which zones were signficant against other zones
  wilcox.result[[i]] <- wilcox$p.value #save pairwise results to a list
}

colnames(kruskal.result) <- c("chi-sq statistic", "df", "pvalue", "test", "variables") 

kruskal.result2 <- kruskal.result %>%
  mutate(significant = ifelse(pvalue < Bonf.alpha, "Yes", "no")) #evaluate if kruskal results signficiant with Bonferroni-corrected p-value

rownames(kruskal.result2) <- variables #add new rownames

names(wilcox.result) <- colnames(data3)

#Save wilcox rank sum results to an excel document
# write.xlsx(wilcox.result, file = "Supp_WilcoxResults.xlsx", row.names = TRUE)

# Which are significant?
kruskal.result2 %>% 
  filter(significant == "Yes") %>%
  rownames()
```

```
## [1] "SCTLD"     "pro"       "syn"       "peuk"      "hbact"     "npoc"     
## [7] "silicate"  "CellRatio"
```

```r
# Graph the significant variables 
dataSig <- data3 %>%
  dplyr::select(date, localtime, station, zonenumber, collab_siteID, FL_Region, FL_region_zone, SCTLD, pro, syn, peuk, hbact, npoc, silicate, CellRatio) %>%
  gather(key = "envdata", value = "value", SCTLD:CellRatio) %>%
  mutate(envdata = factor(envdata, levels = c("SCTLD", "CellRatio", "npoc", "hbact", "pro", "syn", "peuk", "silicate")))

ggplot(dataSig, aes(x = zonenumber, y = value)) + 
  stat_boxplot(geom = "errorbar", width = 0.2) + 
  geom_boxplot(data = dataSig, aes(fill = FL_region_zone), outlier.shape = NA) + 
  geom_point(position = position_jitter(0.2)) +
  labs(title = "Figure 2", x = "FL Coral Reef Zones", fill = "Zone") + 
  scale_x_discrete(limits = rev(levels(dataSig$zonenumber))) +
  scale_fill_manual(values = c("#DC050C", "#FF7F00", "#FDBF6F", "#4EB265", "#1F78B4", "#A6CEE3", "#882E72", "#DE77AE")) +
  theme(panel.grid.major = element_blank(), legend.position = "none") +
  facet_wrap( ~ envdata, scales = "free_y", ncol = 4) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
```

<img src="figures/fig-Environmental data-1.png" width="672" />
