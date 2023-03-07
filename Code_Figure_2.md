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

<img src="figures/fig-2-Environmental data-1.png" width="672" />

# References - R Packages used

```r
knitr::write_bib()
```

```
## Warning in utils::citation(..., lib.loc = lib.loc): no date field in DESCRIPTION
## file of package 'tidyr'
```

```
## @Manual{R-base,
##   title = {R: A Language and Environment for Statistical Computing},
##   author = {{R Core Team}},
##   organization = {R Foundation for Statistical Computing},
##   address = {Vienna, Austria},
##   year = {2022},
##   url = {https://www.R-project.org/},
## }
## 
## @Manual{R-dplyr,
##   title = {dplyr: A Grammar of Data Manipulation},
##   author = {Hadley Wickham and Romain François and Lionel Henry and Kirill Müller and Davis Vaughan},
##   year = {2023},
##   note = {R package version 1.1.0},
##   url = {https://CRAN.R-project.org/package=dplyr},
## }
## 
## @Manual{R-ggplot2,
##   title = {ggplot2: Create Elegant Data Visualisations Using the Grammar of Graphics},
##   author = {Hadley Wickham and Winston Chang and Lionel Henry and Thomas Lin Pedersen and Kohske Takahashi and Claus Wilke and Kara Woo and Hiroaki Yutani and Dewey Dunnington},
##   year = {2022},
##   note = {R package version 3.4.0},
##   url = {https://CRAN.R-project.org/package=ggplot2},
## }
## 
## @Manual{R-tidyr,
##   title = {tidyr: Tidy Messy Data},
##   author = {Hadley Wickham and Davis Vaughan and Maximilian Girlich},
##   year = {2023},
##   note = {https://tidyr.tidyverse.org},
## }
## 
## @Book{ggplot22016,
##   author = {Hadley Wickham},
##   title = {ggplot2: Elegant Graphics for Data Analysis},
##   publisher = {Springer-Verlag New York},
##   year = {2016},
##   isbn = {978-3-319-24277-4},
##   url = {https://ggplot2.tidyverse.org},
## }
```

```r
# what tweak=TRUE does
str(knitr:::.tweak.bib)
```

```
## List of 43
##  $ ade4         : Named chr "  author = {Stéphane Dray and Anne-Béatrice Dufour and Jean Thioulouse and Thibaut Jombart and Sandrine Pavoine"| __truncated__
##   ..- attr(*, "names")= chr "author"
##  $ akima        : Named chr "  author = {H. Akima and Albrecht Gebhardt and Thomas Petzoldt and Martin Maechler},"
##   ..- attr(*, "names")= chr "author"
##  $ ash          : Named chr "  author = {David W. Scott and Albrecht Gebhardt and Stephen Kaluzny},"
##   ..- attr(*, "names")= chr "author"
##  $ bcpa         : Named chr "  author = {Jose Claudio Faria and Clarice Garcia Borges Demetrio},"
##   ..- attr(*, "names")= chr "author"
##  $ BiplotGUI    : Named chr "  author = {Anthony la Grange and  N. J. le Roux and P.J. Rousseeuw and I. Ruts and J. W. Tukey},"
##   ..- attr(*, "names")= chr "author"
##  $ bitops       : Named chr "  author = {Steve Dutky and Martin Maechler},"
##   ..- attr(*, "names")= chr "author"
##  $ cacheSweave  : Named chr "  author = {Roger D. Peng},"
##   ..- attr(*, "names")= chr "author"
##  $ cat          : Named chr "  author = {Ted Harding and Fernando Tusell and Joseph L. Schafer},"
##   ..- attr(*, "names")= chr "author"
##  $ CircStats    : Named chr "  author = {Ulric Lund and Claudio Agostinelli},"
##   ..- attr(*, "names")= chr "author"
##  $ contrast     : Named chr "  author = {Max Kuhn and Steve Weston and Jed Wing and James Forester},"
##   ..- attr(*, "names")= chr "author"
##  $ date         : Named chr "  author = {Terry Therneau and Thomas Lumley and Kjetil Halvorsen and Kurt Hornik},"
##   ..- attr(*, "names")= chr "author"
##  $ digest       : Named chr "  author = {Dirk Eddelbuettel},"
##   ..- attr(*, "names")= chr "author"
##  $ ElemStatLearn: Named chr "  author = {Kjetil Halvorsen},"
##   ..- attr(*, "names")= chr "author"
##  $ epiR         : Named chr "  author = {Mark Stevenson and Telmo Nunes and Cord Heuer and Jonathon Marshall and Javier Sanchez and Ron Thor"| __truncated__
##   ..- attr(*, "names")= chr "author"
##  $ Fahrmeir     : Named chr "  author = {Kjetil Halvorsen},"
##   ..- attr(*, "names")= chr "author"
##  $ flashClust   : Named chr "  author = {Fionn Murtagh and {R development team} and Peter Langfelder},"
##   ..- attr(*, "names")= chr "author"
##  $ foreach      : Named chr "  author = {{Revolution Analytics} and Steve Weston}},"
##   ..- attr(*, "names")= chr "author"
##  $ fortunes     : Named chr "  author = {Achim Zeileis and the R community},"
##   ..- attr(*, "names")= chr "author"
##  $ gee          : Named chr "  author = {Vincent J Carey and Thomas Lumley and Brian Ripley},"
##   ..- attr(*, "names")= chr "author"
##  $ gmodels      : Named chr "  author = {Gregory R. Warnes andBen Bolker and Thomas Lumley and Randall C Johnson and Randall C. Johnson},"
##   ..- attr(*, "names")= chr "author"
##  $ gWidgets     : Named chr "  author = {John Verzani},"
##   ..- attr(*, "names")= chr "author"
##  $ hexbin       : Named chr "  author = {Dan Carr and Nicholas Lewin-Koh and Martin Maechler},"
##   ..- attr(*, "names")= chr "author"
##  $ Hmisc        : Named chr "  author = {Harrell, Jr., Frank E},"
##   ..- attr(*, "names")= chr "author"
##  $ Hmisc        : Named chr "  author = {Frank E. {Harrell, Jr.}},"
##   ..- attr(*, "names")= chr "author"
##  $ leaps        : Named chr "  author = {Thomas Lumley},"
##   ..- attr(*, "names")= chr "author"
##  $ mapproj      : Named chr "  author = {Doug McIlroy and Ray Brownrigg and Thomas P Minka and Roger Bivand},"
##   ..- attr(*, "names")= chr "author"
##  $ maps         : Named chr "  author = {Ray Brownrigg},"
##   ..- attr(*, "names")= chr "author"
##  $ mathgraph    : Named chr "  author = {Patrick J. Burns and Nick Efthymiou and Claus Dethlefsen},"
##   ..- attr(*, "names")= chr "author"
##  $ oz           : Named chr "  author = {Bill Venables and Kurt Hornik},"
##   ..- attr(*, "names")= chr "author"
##  $ pbivnorm     : Named chr "  author = {Alan Genz and Brenton Kenkel},"
##   ..- attr(*, "names")= chr "author"
##  $ pscl         : Named chr "  author = {Simon Jackman and Alex Tahk and Achim Zeileis and Christina Maimone and Jim Fearon},"
##   ..- attr(*, "names")= chr "author"
##  $ quadprog     : Named chr "  author = {Berwin A. Turlach and Andreas Weingessel},"
##   ..- attr(*, "names")= chr "author"
##  $ R2SWF        : Named chr "  author = {Yixuan Qiu and Yihui Xie and Cameron Bracken},"
##   ..- attr(*, "names")= chr "author"
##  $ R2WinBUGS    : Named chr "  author = {Andrew Gelman and Sibylle Sturtz and Uwe Ligges and Gregor Gorjanc and Jouni Kerman},"
##   ..- attr(*, "names")= chr "author"
##  $ randomForest : Named chr "  author = {Leo Breiman and Adele Cutler and Andy Liaw and Matthew Wiener},"
##   ..- attr(*, "names")= chr "author"
##  $ rgl          : Named chr "  author = {Daniel Adler and Duncan Murdoch},"
##   ..- attr(*, "names")= chr "author"
##  $ RgoogleMaps  : Named chr "  author = {Markus Loecher},"
##   ..- attr(*, "names")= chr "author"
##  $ rms          : Named chr "  author = {Frank E. {Harrell, Jr.}},"
##   ..- attr(*, "names")= chr "author"
##  $ robustbase   : Named chr "  author = {Valentin Todorov and Andreas Ruckstuhl and Matias Salibian-Barrera and Tobias Verbeke and Manuel Ko"| __truncated__
##   ..- attr(*, "names")= chr "author"
##  $ RODBC        : Named chr "  author = {Brian Ripley and Michael Lapsley},"
##   ..- attr(*, "names")= chr "author"
##  $ Sleuth2      : Named chr "  author = {F. L. Ramsey and D. W. Schafer and Jeannie Sifneos and Berwin A. Turlach},"
##   ..- attr(*, "names")= chr "author"
##  $ sm           : Named chr "  author = {Adrian Bowman and Adelchi Azzalini},"
##   ..- attr(*, "names")= chr "author"
##  $ tuneR        : Named chr "  author = {Uwe Ligges},"
##   ..- attr(*, "names")= chr "author"
```
