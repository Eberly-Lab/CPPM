setwd("Set/Path")


library(ggplot2)
library(ggpubr)
library(microeco)
library(file2meco)
library(mecodev)
library(dplyr)
library(tidyr)
library(vegan)
library(lavaan)  
library(semPlot) 
library(RColorBrewer)
theme_set(theme_classic())


# The following code is for Fig. 2 and Fig. S1
m16S <- readRDS("m16S_combined.RDS") 

m16S$cal_betadiv(unifrac = TRUE)

sd <- as.data.frame(m16S$sample_table)
a_div <- trans_alpha$new(dataset = m16S, group = "Year", by_group = "Location")
a_df <- as.data.frame(a_div$data_alpha)
a_df <- a_df[, c(1:3,11:29,37)]

a_dfw <- pivot_wider(a_df, names_from = "Measure", values_from = "Value")

# Recode categorical as numeric
a_dfwn <- a_dfw
a_dfwn$Seed_num <- recode(a_dfwn$Seeding_Rate, "High" = 1, "Low" = 0)
a_dfwn$N_num <- recode(a_dfwn$Nitrogen, "High" = 1, "Low" = 0)
a_dfwn$Fun_num <- recode(a_dfwn$Fungicide, "Fungicide" = 1, "No Fungicide" = 0)
a_dfwn$Weeds_num <- recode(a_dfwn$B.tectorum, "B. tectorum" = 1, "No B. tectorum" = 0)
a_dfwn$F.p <- a_dfwn$F.pseudograminearum_count
a_dfwn$Observed_div <- a_dfwn$Observed/100
a_dfwn <- as.data.frame(a_dfwn)
rownames(a_dfwn) <- a_dfwn$Sample

otu_df <- as.data.frame(m16S$otu_table)
bray_dist <- vegdist(t(otu_df), method = "bray")

# Perform PCoA on beta diversity distances
pcoa_bray <- cmdscale(bray_dist, k = 5, eig = TRUE)
pcoa_bray$points[, 1]
m16S$beta_diversity$unwei_unifrac[,1]

# Extract PCoA axes (these will be our dependent variables in SEM)
pcoa_axes <- data.frame(
  Sample = sd$Sample,
  PC1 = pcoa_bray$points[, 1],
  PC2 = pcoa_bray$points[, 2],
  wuni = m16S$beta_diversity$wei_unifrac[, 1],
  uni = m16S$beta_diversity$unwei_unifrac[, 1]
  )

# Combined dataframe with sample data, alpha diversity, and PCoA axes
am_df <- merge(a_dfwn, pcoa_axes, by = 0)

# Add network data
bz22f_hub <- read.csv("Networks/Bozeman 2022 Fusarium_Degree_top20 default node.csv")
bz22nf_hub <- read.csv("Networks/Bozeman 2022 No Fusarium_Degree_top20 default node.csv")
bz23f_hub <- read.csv("Networks/Bozeman 2023 Fusarium_Degree_top20 default node.csv")
bz23nf_hub <- read.csv("Networks/Bozeman 2023 No Fusarium_Degree_top20 default node.csv")
sarc22f_hub <- read.csv("Networks/Huntley 2022 Fusarium_Degree_top20 default node.csv")
sarc22nf_hub <- read.csv("Networks/Huntley 2022 No Fusarium_Degree_top20 default node.csv")
sarc23f_hub <- read.csv("Networks/Huntley 2023 Fusarium_Degree_top20 default node.csv")
sarc23nf_hub <- read.csv("Networks/Huntley 2023 Fusarium_Degree_top20 default node.csv")

bz22f_hub$Location <- "BZ22f"
bz22nf_hub$Location <- "BZ22nf"
bz23f_hub$Location <- "BZ23f"
bz23nf_hub$Location <- "BZ23nf"
sarc22f_hub$Location <- "SARC22f"
sarc22nf_hub$Location <- "SARC22nf"
sarc23f_hub$Location <- "SARC23f"
sarc23nf_hub$Location <- "SARC23nf"

fus_hub <- rbind(bz22f_hub, bz23f_hub, sarc22f_hub, sarc23f_hub)
nfus_hub <- rbind(bz22nf_hub, bz23nf_hub, sarc22nf_hub, sarc23nf_hub)

hub <- rbind(bz22f_hub, bz23f_hub, sarc22f_hub, sarc23f_hub, bz22nf_hub, bz23nf_hub, sarc22nf_hub, sarc23nf_hub)

fus_asv <- fus_hub$id
nfus_asv <- nfus_hub$id
asv <- hub$id

m16Sn <- trans_norm$new(dataset = m16S)
m16Sn <- m16Sn$norm(method = "clr")

fus_tk <- rownames(m16Sn$otu_table) %in% fus_asv
nfus_tk <- rownames(m16Sn$otu_table) %in% nfus_asv
hub_tk <- rownames(m16Sn$otu_table) %in% asv

fus_otu <- m16Sn$otu_table
fus_otu <- fus_otu[fus_tk,]
fus_otu <- t(fus_otu)
fus_otu <- data.frame(fus_otu, check.names = FALSE)
fus_otu$F_Hub <- rowSums(fus_otu)
fus_otu <- fus_otu[,c(79,1:78)]

nfus_otu <- m16Sn$otu_table
nfus_otu <- nfus_otu[nfus_tk,]
nfus_otu <- t(nfus_otu)
nfus_otu <- data.frame(nfus_otu, check.names = FALSE)
nfus_otu$NF_Hub <- rowSums(nfus_otu)
nfus_otu <- nfus_otu[,c(80,1:79)]

otu <- m16Sn$otu_table
otu <- otu[hub_tk,]
otu <- t(otu)
otu <- data.frame(otu, check.names = FALSE)
otu$Hub <- (rowSums(otu))
otu <- otu[,c(134,1:133)]

rownames(am_df) <- am_df$Row.names
am_df[,1] <- NULL

am_merge <- merge(am_df, fus_otu[,c("F_Hub", "545864")], by = 0, all.x = TRUE)
rownames(am_merge) <- am_merge$Row.names
am_merge[,1] <- NULL
am_merge <- merge(am_merge, nfus_otu[,c("NF_Hub", "545864")], by = 0, all.x = TRUE)
am_merge <- am_merge[,-c(2,43,45)]
rownames(am_merge) <- am_merge$Row.names
am_merge[,1] <- NULL
am_merge <- merge(am_merge, otu[,c("Hub", "545864")], by = 0, all.x = TRUE)
rownames(am_merge) <- am_merge$Row.names
am_merge <- am_merge[,-c(1,44)]

# Full model
l.mod <- ' Observed_div ~ Loc_Yr + F.p + B.tectorum + Seeding_Rate + nitrate + phosphorus + OM
PC1 ~ F.p + Loc_Yr + Seeding_Rate + B.tectorum + Fungicide + nitrate + phosphorus + OM
uni ~ F.p + Loc_Yr + Seeding_Rate + B.tectorum + Fungicide + nitrate + phosphorus + OM
Hub ~ F.p + Loc_Yr + Seeding_Rate + B.tectorum + Fungicide + nitrate + phosphorus + OM
'
# Reduced model
l.mod1 <- ' Observed_div ~ Loc_Yr + F.p
PC1 ~ Loc_Yr + F.p
uni ~ Loc_Yr + OM + F.p
Hub ~ F.p
'

sem1 <- sem(l.mod, data = am_merge)
sem1c <- sem(l.mod, data = am_merge, group.equal = c("intercepts", "regressions"))
sem2 <- sem(l.mod1, data = am_merge)
anova(sem1, sem2)

# Use model 1 for the complete model and Reduced model for plot
summary(sem1, standardize = T, rsq = T, fit.measures=T)
summary(sem2, standardize = T, rsq = T, fit.measures=T)

# Capture full summary output
sink("sem_model_summary.txt")
summary(sem1, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
sink()

# Capture reduced model summary output
sink("sem_reduced_model_summary.txt")
summary(sem2, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
sink()

param_df <- parameterEstimates(sem1, standardized = TRUE)
write.csv(param_df, "SEM parameters full model.csv")

param_df1 <- parameterEstimates(sem2, standardized = TRUE)
write.csv(param_df1, "SEM parameters reduced model.csv")

fitMeasures(sem2)
varTable(sem1)
anova(sem1, sem2)

semPaths(sem1, what = "std", layout = "spring", style = "OpenMx", edge.label.cex = 0.8, curvePivot = TRUE, residuals = FALSE,
         groups = list("Soil Properties" = c("ntr", "phs", "OM"), "Management Practices" = c("S_R",
            "Fng", "Ntr"), "Bacterial Community" = c("Ob_", "PC1", "uni", "Hub"), "Treatment" = 
              c("Fsr", "Wds", "L_Y")),  inheritColor = TRUE, pastel = TRUE, shapeMan = "circle", cut = 0.1)

pdf("Figures/FigS1.pdf", width = 10, height = 8)
semPaths(sem1, what = "std", layout = "spring", style = "OpenMx", edge.label.cex = 0.8, curvePivot = TRUE, residuals = FALSE,
         groups = "manifest", pastel = TRUE, shapeMan = "circle", cut = 0.1, edge.label.position = 0.6)
dev.off()


# Add p-values
table1<-parameterEstimates(sem1,standardized=TRUE) #  %>%  head(49)
table1 <- subset(table1, pvalue <=0.05)

b1<-gettextf('%.3f \n p=%.3f', table1$std.all, digits=table1$pvalue)


# Reduced model plot
semPaths(sem2, what = "std", layout = "spring", style = "OpenMx", edge.label.cex = 0.8, curvePivot = TRUE, residuals = FALSE,
         groups = "manifest", pastel = TRUE, shapeMan = "circle", cut = 0.1)

# Add p-values
table2<-parameterEstimates(sem2,standardized=TRUE) 
table2[is.na(table2)] <- 1

# Remove self comparisons, i.e., obs ~~ obs
table2 <- table2[-c(9:12,19,22,24),]

b2<-gettextf('%.3f \n p=%.3f', table2$std.all, digits=table2$pvalue)



# Parameter estimates plot
p1 <- parameterEstimates(sem1, standardized = TRUE) %>%
  filter(op == "~") %>%
  ggplot(aes(x = rhs, y = std.all, fill = lhs)) +
  geom_col(position = "dodge") +
  facet_wrap(~lhs, scales = "free") +
  theme_minimal() +
  labs(x = "Predictor", y = "Standardized Path Coefficient") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
p1

pdf("Figures/SEM Parameters Full Model.pdf", width = 10, height = 8)
p1
dev.off()

p2 <- parameterEstimates(sem2, standardized = TRUE) %>%
  filter(op == "~") %>%
  ggplot(aes(x = rhs, y = std.all, fill = lhs)) +
  geom_col(position = "dodge") +
  facet_wrap(~lhs, scales = "free") +
  theme_minimal() +
  labs(x = "Predictor", y = "Standardized Path Coefficient") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
p2

# Fig. 2
library(grid)
library(ggplotify)
library(ggpubr)

grec <- recordPlot(semPaths(sem2, what = "std", edgeLabels = b2, style="OpenMx",edge.label.cex = 0.65, edge.label.position = 0.72, curvePivot = TRUE, layout = 'tree',
                            intercepts=FALSE, residuals=FALSE, groups = "manifest", inheritColor = TRUE, pastel = TRUE, shapeMan = "circle", 
                            cut = 0.05, alpha = 0.01, fade = TRUE))

ggarrange(grec, p1, nrow = 2, ncol = 1, widths = c(1,0.5))

pdf("Figures/Fig2.pdf", width = 10, height = 10)
ggarrange(grec, p1, nrow = 2, ncol = 1, heights = c(1,0.75), labels = c("a", "b"))
dev.off()

