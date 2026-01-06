setwd("your/path/here")


library(phyloseq)
library(dplyr)
library(vegan)
library(ggplot2)
library(ggpubr)
library(magrittr)
library(microbiomeutilities)
library(microbiome)
library(car)
library(plyr)
library(paletteer)

# Import phyloseq objects generated in CPPM Emu.R
ps_bz_22 <- readRDS("emu_phyloseq_bz_2022.RDS")
ps_bz_23 <- readRDS("emu_phyloseq_bz_2023.RDS")
ps_sarc_22 <- readRDS("emu_phyloseq_SARC_2022.RDS")
ps_sarc_23 <- readRDS("emu_phyloseq_SARC_2023.RDS")


# Check number of common and unique sequences
bz22n <- taxa_names(ps_bz_22)
bz23n <- taxa_names(ps_bz_23)
sarc22n <- taxa_names(ps_sarc_22)
sarc23n <- taxa_names(ps_sarc_23)

v1 <- list("BZ22" = bz22n, "BZ23" = bz23n, "SARC22" = sarc22n, "SARC23" = sarc23n)

library(ggVennDiagram)
ggVennDiagram(v1, label_percent_digit = 2)

# Make combined phyloseq object
otu_bz_22 <- otu_table(ps_bz_22)
otu_bz_23 <- otu_table(ps_bz_23)
otu_sarc_22 <- otu_table(ps_sarc_22)
otu_sarc_23 <- otu_table(ps_sarc_23)

tax_bz_22 <- tax_table(ps_bz_22)
tax_bz_23 <- tax_table(ps_bz_23)
tax_sarc_22 <- tax_table(ps_sarc_22)
tax_sarc_23 <- tax_table(ps_sarc_23)

sd_bz_22 <- sample_data(ps_bz_22)
sd_bz_23 <- sample_data(ps_bz_23)
sd_sarc_22 <- sample_data(ps_sarc_22)
sd_sarc_23 <- sample_data(ps_sarc_23)

# Add Year and Location
sd_bz_22$Year <- "2022"
sd_bz_23$Year <- "2023"
sd_sarc_22$Year <- "2022"
sd_bz_22$Location <- "Bozeman"
sd_bz_23$Location <- "Bozeman"
sd_sarc_22$Location <- "Huntley"
sd_sarc_23$Location <- "Huntley"

# Append row names in sd and column names in tax_tab and otu_tab to retain unique taxa
sd_BZ_22rn <- rownames(sd_bz_22)
sd_BZ_23rn <- rownames(sd_bz_23)
sd_SARC_22rn <- rownames(sd_sarc_22)
sd_SARC_23rn <- rownames(sd_sarc_23)

sd_BZ_22rn <- paste0(sd_BZ_22rn, "_BZ_22")
sd_BZ_23rn <- paste0(sd_BZ_23rn, "_BZ_23")
sd_SARC_22rn <- paste0(sd_SARC_22rn, "_SARC_22")
sd_SARC_23rn <- paste0(sd_SARC_23rn, "_SARC_23")

rownames(sd_bz_22) <- sd_BZ_22rn
rownames(sd_bz_23) <- sd_BZ_23rn
rownames(sd_sarc_22) <- sd_SARC_22rn
rownames(sd_sarc_23) <- sd_SARC_23rn

otu_bz_22_cn <- colnames(otu_bz_22)
otu_bz_23_cn <- colnames(otu_bz_23)
otu_sarc_22_cn <- colnames(otu_sarc_22)
otu_sarc_23_cn <- colnames(otu_sarc_23)

otu_bz_22_cn <- paste0(otu_bz_22_cn, "_BZ_22")
otu_bz_23_cn <- paste0(otu_bz_23_cn, "_BZ_23")
otu_sarc_22_cn <- paste0(otu_sarc_22_cn, "_SARC_22")
otu_sarc_23_cn <- paste0(otu_sarc_23_cn, "_SARC_23")

colnames(otu_bz_22) <- otu_bz_22_cn
colnames(otu_bz_23) <- otu_bz_23_cn
colnames(otu_sarc_22) <- otu_sarc_22_cn
colnames(otu_sarc_23) <- otu_sarc_23_cn


ps1 <- phyloseq(otu_bz_22, tax_bz_22, sd_bz_22)
ps2 <- phyloseq(otu_bz_23, tax_bz_23, sd_bz_23)
ps3 <- phyloseq(otu_sarc_22, tax_sarc_22, sd_sarc_22)
ps4 <- phyloseq(otu_sarc_23, tax_sarc_23, sd_sarc_23)

ps <- merge_phyloseq(ps1, ps2, ps3, ps4)
ps

table(tax_table(ps)[, "Phylum"], exclude = NULL)
ps1 <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
table(tax_table(ps1)[, "Phylum"], exclude = NULL)

saveRDS(ps1, "merged_phyloseq.RDS")
sd <- sample_data(ps)
sd <- data.frame(sd)
write.csv(sd, "combined_sample_data.csv")

ps <- phyloseq(otu_table(ps), tax_table(ps), sample_data(sd), refseq(ps))
ps

sample_data(ps)$Loc_Yr <- paste0(sample_data(ps)$Location, sep="_", sample_data(ps)$Year)
sample_data(ps)

saveRDS(ps, "merged_phyloseq.RDS")

ps_ref <- add_refseq(ps, tag = "ASV")
ps_ref

##########################################
# Start Here With Merged Phyloseq Object #
##########################################

library(microeco)
library(file2meco)
library(mecodev)
theme_set(theme_classic())

ps <- readRDS("merged_phyloseq.RDS")

m16S <- phyloseq2meco(ps)
saveRDS(m16S, "m16S_combined.RDS")
m16S <- readRDS("m16S_combined.RDS")
sd <- m16S$sample_table


# Rarefaction curve
r1 <- trans_rarefy$new(m16S, alphadiv = "Observed", depth = c(0, 10, 50, 500, 2000, 4000, 6000, 8000, 10000, 20000))
p1 <- r1$plot_rarefy(color = "Location", show_point = FALSE, add_fitting = FALSE)
p1
p1 + theme(legend.position="none") 
saveRDS(r1, "rarefy_16S.RDS")

#m16S$rarefy_samples(sample.size = 5000)

# Permanova
m16S$cal_betadiv()
m16S$tidy_dataset()

p_mod <- trans_beta$new(dataset = m16S, measure = "bray")
p_mod$cal_manova(manova_set = "Location*Year+Weeds*Fusarium+Nitrogen+Fungicide+Seeding_Rate+Cochliob", permutations = 9999) 
p_mod$res_manova
# Significant Location X Year interaction
# This is justification to evaluate treatment effects by individual location-year

mod_df <- as.data.frame(p_mod$res_manova)
write.csv(mod_df, "Permanova.csv")

ggboxplot(sd, x = "Weeds", y = "Fusarium") 

# by location year
bz22_16S <- clone(m16S)
bz23_16S <- clone(m16S)
sarc22_16S <- clone(m16S)
sarc23_16S <- clone(m16S)

bz22_16S$sample_table <- subset(bz22_16S$sample_table, Loc_Yr == "Bozeman_2022")
bz23_16S$sample_table <- subset(bz23_16S$sample_table, Loc_Yr == "Bozeman_2023")
sarc22_16S$sample_table <- subset(sarc22_16S$sample_table, Loc_Yr == "Huntley_2022")
sarc23_16S$sample_table <- subset(sarc23_16S$sample_table, Loc_Yr == "Huntley_2023")
bz22_16S$tidy_dataset()
bz23_16S$tidy_dataset()
sarc22_16S$tidy_dataset()
sarc22_16S$tidy_dataset()

p_bz22_16S <- trans_beta$new(dataset = bz22_16S, measure = "bray")
p_bz23_16S <- trans_beta$new(dataset = bz23_16S, measure = "bray")
p_sarc22_16S <- trans_beta$new(dataset = sarc22_16S, measure = "bray")
p_sarc23_16S <- trans_beta$new(dataset = sarc23_16S, measure = "bray")

p_bz22_16S$cal_manova(manova_set = "Weeds*Fusarium*Nitrogen*Fungicide+Seeding_Rate+Cochliob", permutations = 9999) 
p_bz23_16S$cal_manova(manova_set = "Weeds*Fusarium*Nitrogen*Fungicide+Seeding_Rate+Cochliob", permutations = 9999) 
p_sarc22_16S$cal_manova(manova_set = "Weeds*Fusarium*Nitrogen*Fungicide+Seeding_Rate+Cochliob", permutations = 9999) 
p_sarc23_16S$cal_manova(manova_set = "Weeds*Fusarium*Nitrogen*Fungicide+Seeding_Rate+Cochliob", permutations = 9999) 

p_bz22_16S$res_manova
p_bz23_16S$res_manova
p_sarc22_16S$res_manova
p_sarc23_16S$res_manova


# Alpha Diversity
a_div <- trans_alpha$new(dataset = m16S, group = "Year", by_group = "Location")
a_div$cal_diff(method = "KW")

a_div$plot_alpha(measure = "InvSimpson")
a_div$plot_alpha(measure = "Shannon")
a_div$plot_alpha(measure = "Coverage")

# Linear model
lm_df <- as.data.frame(m16S$alpha_diversity)
sd <- as.data.frame(m16S$sample_table)
lm_df <- merge(sd, lm_df, by = 0)

library(lme4)
library(car)
ex_mod <- lmer(Chao1 ~ Location*Year+Seeding_Rate+Weeds+Nitrogen+Fungicide+Fusarium + (1|Block), data = lm_df)
summary(ex_mod)
Anova(ex_mod, type="II")

mod_df <- data.frame(Anova(ex_mod, type="II"))
mod_df$Factor <- rownames(mod_df)
mod_df <- mod_df[,c(4,1,2,3)]
write.csv(mod_df, "alph div model.csv")



# Differential abundance analysis
t1 <- trans_diff$new(dataset = m16S, method = "ALDEx2_kw", group = "Weeds", taxa_level = "Genus")
t1 <- trans_diff$new(dataset = m16S, method = "ALDEx2_kw", group = "FusariumPresence", taxa_level = "Family")
# see t1$res_diff for the result
range(t1$res_diff$P.adj)
t1$plot_diff_bar(use_number = 1:20, width = 0.8)
t1$plot_diff_abund()

dim(m16S$otu_table)

m16_bz22 <- phyloseq2meco(bz22)
m16_bz23 <- phyloseq2meco(bz23)
m16_sarc22 <- phyloseq2meco(sarc22)
m16_sarc23 <- phyloseq2meco(sarc23)

ald_bz22 <- trans_diff$new(dataset = m16_bz22, method = "ALDEx2_kw", group = "FusariumPresence", taxa_level = "Genus")
ald_bz23 <- trans_diff$new(dataset = m16_bz23, method = "ALDEx2_kw", group = "FusariumPresence", taxa_level = "Genus")
ald_sarc22 <- trans_diff$new(dataset = m16_sarc22, method = "ALDEx2_kw", group = "FusariumPresence", taxa_level = "Family")
ald_sarc23 <- trans_diff$new(dataset = m16_sarc23, method = "ALDEx2_kw", group = "FusariumPresence", taxa_level = "Family")

range(ald_bz22$res_diff$P.adj)
ald_bz22$plot_diff_abund()
ald_bz23$plot_diff_abund()
ald_sarc22$plot_diff_abund()
ald_sarc23$plot_diff_abund()

# SEM of alpha diversity metrics
library(piecewiseSEM)
library(tidyr)
library(nlme)

a_df <- as.data.frame(a_div$data_alpha)
a_df <- a_df[, c(2,3,12:24,36)]

a_dfw <- pivot_wider(a_df, names_from = "Measure", values_from = "Value")

sem_mod <- psem(
  lm(Observed ~ Loc_Yr + Seeding_Rate + Nitrogen + Weeds + Fungicide + FusariumPresence, data = a_dfw),
  lm(Chao1 ~ Loc_Yr + Seeding_Rate + Nitrogen + Weeds + Fungicide + FusariumPresence, data = a_dfw),
  lm(Shannon ~ Loc_Yr + Seeding_Rate + Nitrogen + Weeds + Fungicide + FusariumPresence, data = a_dfw),
  lm(Fisher ~ Loc_Yr + Seeding_Rate + Nitrogen + Weeds + Fungicide + FusariumPresence, data = a_dfw)
)

summary(sem_mod, .progressBar = TRUE)

sem_mod1 <- psem(
  lme(Observed ~ Seeding_Rate + Nitrogen + Weeds + Fungicide + FusariumPresence, data = a_dfw),
  lme(Chao1 ~ Seeding_Rate + Nitrogen + Weeds + Fungicide + FusariumPresence, data = a_dfw),
  lme(Shannon ~ Seeding_Rate + Nitrogen + Weeds + Fungicide + FusariumPresence, data = a_dfw),
  lme(Fisher ~ Seeding_Rate + Nitrogen + Weeds + Fungicide + FusariumPresence, data = a_dfw)
)
summary(sem_mod1, .progressBar = TRUE)

library(lavaan)
library(dplyr)
# Categorical as numeric
a_dfwn <- a_dfw
a_dfwn$Seed_num <- dplyr::recode(a_dfwn$Seeding_Rate, "High" = 1, "Low" = 0)
a_dfwn$N_num <- dplyr::recode(a_dfwn$Nitrogen, "High" = 1, "Low" = 0)
a_dfwn$Fun_num <- dplyr::recode(a_dfwn$Fungicide, "Fungicide" = 1, "No Fungicide" = 0)
a_dfwn$Weeds_num <- dplyr::recode(a_dfwn$Weeds, "Brome" = 1, "No Brome" = 0)
a_dfwn$Obs_div <- a_dfwn$Observed/100

l.mod <- ' Observed ~ Seeding_Rate + Nitrogen + Weeds + Fungicide + FusariumPresence
Chao1 ~  Seeding_Rate + Nitrogen + Weeds + Fungicide + FusariumPresence
Shannon ~  Seeding_Rate + Nitrogen + Weeds + Fungicide + FusariumPresence
Fisher ~  Seeding_Rate + Nitrogen + Weeds + Fungicide + FusariumPresence
'
sem_mod1 <- cfa(l.mod, data = a_dfwn, group = "Loc_Yr")


l.mod <- ' Observed ~ Loc_Yr + Seed_num + N_num + Weeds_num + Fun_num + Fusarium
Fusarium ~ Loc_Yr + Seed_num + N_num + Weeds_num + Fun_num + Cochliob
'

l.mod <- ' Obs_div~ Fusarium
Fusarium ~ Loc_Yr + Seed_num + N_num + Weeds_num + Fun_num
'

l.mod1 <- ' Obs_div ~ Loc_Yr + Fusarium
Fusarium ~ Loc_Yr + Seeding_Rate + Nitrogen + Weeds + Fungicide

'

sem1 <- sem(l.mod, data = a_dfwn)
sem2 <- sem(l.mod1, data = a_dfwn)

summary(sem2, standardize = T, rsq = T, fit.measures=T)
fitMeasures(sem2)
varTable(sem1)

anova(sem1, sem2)

# Plots
library(MicrobiotaProcess)

sdf <- data.frame(sample_data(ps))
sdf$Loc_Yr <- paste(sdf$Location, sdf$Year, sep="_")
sample_data(ps) <- sample_data(sdf)

mpse1 <- ps %>% as.MPSE()

mp.rar <- mpse1 %>% mp_rrarefy()

mp.rar %<>% 
  mp_cal_alpha()

a_df <- as.data.frame(a_div$data_alpha)

cov_df <- subset(a_df, Measure == "Observe" | Measure == "InvSimpson")
cov_df <- cov_df[, c(1:3)]
cov_df <- pivot_wider(cov_df, names_from = "Measure", values_from = "Value")
cov_df <- as.data.frame(cov_df)
rownames(cov_df) <- cov_df$Sample  
cov_df <- cov_df[,c(2,3)]

mp.rar$Coverage <- cov_df$Coverage
mp.rar$InvSimpson <- cov_df$InvSimpson

mp.rar  <- merge(mp.rar, cov_df, by = "Sample")
mp.rar <- as_tibble(mp.rar)
mp.rar$Coverage

p1 <- mp.rar %>% 
  mp_plot_alpha(
    .group=c(Loc_Yr), 
    .alpha=c(Observe, InvSimpson),
    step_increase = 0.1
  ) + theme(legend.position="none") +
  geom_signif(test="wilcox.test", textsize=4) +
  theme_classic() + theme(text = element_text(size = 14), 
                          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
p1 <- p1 + theme(legend.position="none")

pdf("Figures/Alpha Diversity Loc_Yr.pdf", width = 10, height = 8)
p1
dev.off()


p2 <- mpse1 %>% 
  mp_plot_alpha(
    .group=c(Treatment4, Loc_Yr), 
    .alpha=c(Chao1),
    step_increase = 0.1
  ) + theme(legend.position="none") +
  geom_signif(comparisons=list(c("Bozeman", "Huntley")), test="wilcox.test", textsize=4)
p2


# Differences in alpha and beta diversity were significant between locations and years which is
# justification to evaluate location years separately

#### Make bar plots and ordination plots for multi panel Fig. 1.
# Remaining analysis done on individual location years


ps <- readRDS("merged_phyloseq.RDS")
sd <- sample_data(ps)
sd <- data.frame(sd)
write.csv(sd, "combined_sample_data.csv")

sd1 <- read.csv("combined_sample_data.csv", header = TRUE, row.names = 1)



############################
# Beta Diversity Bar Plots #
############################

m16S$sample_table$Loc_Yr <- paste(m16S$sample_table$Location, m16S$sample_table$Year, sep="_")
bpp <- trans_abund$new(dataset = m16S, taxrank = "Phylum", ntaxa = 10, groupmean = "Loc_Yr") 
bpp <- bpp$plot_bar(others_color = "grey70", legend_text_italic = FALSE) + 
  theme_classic() + theme(text = element_text(size = 14), 
  axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
bpp

pdf("Figures/Beta Diversity Loc Yr.pdf", width = 10, height = 8)
bpp
dev.off()

bpf <- trans_abund$new(dataset = m16S, taxrank = "Family", ntaxa = 29, groupmean = "Loc_Yr") 
bpf1 <- bpf$plot_bar(others_color = "grey70", legend_text_italic = FALSE) + 
  theme_classic() + theme(text = element_text(size = 14), 
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
bpf1 

pdf("Figures/Beta Diversity Family Loc Yr.pdf", width = 10, height = 8)
bpf1 + theme_classic() + theme(text = element_text(size = 14))
dev.off()

bpbf <- trans_abund$new(dataset = m16S, taxrank = "Family", ntaxa = 15) 
bpbf$plot_box(group = "Loc_Yr") 

t1 <- trans_abund$new(dataset = m16S, taxrank = "Class", ntaxa = 20, show = 0, high_level = "Phylum", high_level_fix_nsub = 3, delete_taxonomy_prefix = FALSE)
t1$plot_bar(ggnested = TRUE, high_level_add_other = TRUE, xtext_angle = 30, facet = c("Treatment1", "Loc_Yr"))


###################
# microViz method #
###################
library(microViz)
library(GUniFrac)
library(RColorBrewer)
library(phyloseq)
library(microbiome)
library(ggpubr)

psv <- tax_fix(ps)
#tax_fix_interactive(ps)

sdf <- data.frame(sample_data(psv))
sd$Loc_Yr <- paste(sd$Location, sdf$Year, sep="_")

# For multiple variables
variables_to_convert <- c("nitrate", "phosphorus", "potassium", "OM")

for (var in variables_to_convert) {
  sdf[[var]] <- as.factor(sdf[[var]])
}

sample_data(psv) <- sample_data(sdf)

ord1 <- psv %>%
  tax_filter(min_prevalence = 2.5 / 100, verbose = FALSE) %>%
  
  dist_calc("bray") %>%
  ord_calc(method = "NMDS") %>%
  ord_plot(alpha = 0.6, size = 2, color = "Loc_Yr", auto_caption = NA) +
  scale_color_brewer(name = "Location", palette = "Dark2") +
  theme_classic(12) + stat_ellipse(aes(color = Loc_Yr))
ord1

# RDA
rda_a <- psv %>% tax_fix(unknowns = c("uncultured"))
rda_a <- subset_taxa(rda_a, !is.na(Genus) & !Genus %in% c("uncultured"))

rda1 <- rda_a %>%
  tax_transform("clr", rank = "Phylum") %>%
  ord_calc(
    method = "RDA", # Note: you can specify RDA explicitly, and it is good practice to do so, but microViz can guess automatically that you want an RDA here (helpful if you don't remember the name?)
    scale_cc = FALSE # doesn't make a difference
  ) %>% ord_plot(color = "Loc_Yr", size = 2, alpha = 0.5, auto_caption = NA) + 
  theme_classic() + stat_ellipse(aes(color = Loc_Yr))
rda1

# Multi panel figure
g1 <- ggarrange(p1, ord1, ncol = 2, labels = c("a", "b")) +  
  theme(plot.margin = margin(0.1,0.1,0.5,0.1, "cm")) 
g1

g2 <- ggarrange(bpp, bpf1, ncol = 2, widths = c(1, 1.2), labels = c("c", "d")) +  
  theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm")) 
g2

ggarrange(g1, g2, nrow = 2, widths = c(1,1))

pdf("Figures/Fig1.pdf", width = 12, height = 10)
ggarrange(g1, g2, nrow = 2)
dev.off()



# Upset Plot by Location
brome <- m16S$merge_samples(group = "Weeds")
t1 <- trans_venn$new(brome, ratio = "numratio")
t1$plot_venn()

fung <- m16S$merge_samples(group = "Fungicide")
t2 <- trans_venn$new(fung, ratio = "numratio")
t2$plot_venn()

trt <- m16S$merge_samples(group = "Treatment2")
t3 <- trans_venn$new(trt, ratio = "numratio")
g1 <- t3$plot_bar(left_plot = TRUE, bottom_height = 0.5, left_width = 0.15, up_bar_fill = "grey50", left_bar_fill = "grey50", bottom_point_color = "black")
g1

pdf("Figures/Brome Fungicide Upset Plot.pdf", width = 10, height = 6)
g1
dev.off()

sd$Treatment8 <- paste(sd$Weeds, sd$Fusarium)

# Fig. 3: Correlation plot
m16S <- readRDS("microeco 16S.RDS") 

sub_16S <- clone(m16S)
sub_16S <- sub_16S$filter_taxa(rel_abund = 0.004, freq = 0.5)

sub_16S$sample_table$N_rate <- dplyr::recode(sub_16S$sample_table$Nitrogen, "High" = 1, "Low" = 0)
sub_16S$sample_table$Seeding_Rate <- dplyr::recode(sub_16S$sample_table$Seeding_Rate, "High" = 1, "Low" = 0)
sub_16S$sample_table$B.tectorum <- dplyr::recode(sub_16S$sample_table$B.tectorum, "B. tectorum" = 1, "No B. tectorum" = 0)

sub_16S$sample_table[,c(9,11,13,16,17, 34)]
cor_env <- trans_env$new(dataset = sub_16S, add_data = sub_16S$sample_table[, c(9,11,13,16,17,34)])

cor_env$cal_cor(use_data = "Species")
cor_env$plot_cor()

pdf("Figures/Fig3.pdf", width = 8, height = 10)
cor_env$plot_cor()
dev.off()

cor_tax <- cor_env$res_cor
write.csv(cor_tax, "taxa correlation.csv")
mas <- subset(cor_env$dataset$tax_table, Genus == "g__Massilia")
baci <- subset(m16S$tax_table, Genus == "g__Bacillus")

sub_16S$sample_table[,c(17)]
cor_env1 <- trans_env$new(dataset = sub_16S, add_data = sub_16S$sample_table[, c(17), drop = FALSE])
cor_env1$cal_cor(use_data = "Species")
cor_env1$res_cor <- subset(cor_env1$res_cor, Correlation > 0.10 | Correlation < -0.10)
cor_env1$plot_cor()







