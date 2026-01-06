setwd("set/path")


library(phyloseq)
library(ggplot2)
library(ggpubr)
library(microeco)
library(file2meco)
library(mecodev)
library(meconetcomp)
library(magrittr)
library(microViz)
library(multcompView)

theme_set(theme_classic())

# Compare networks to determine if there is a fusarium brome interaction
m16 <- readRDS("m16S_combined.RDS") # load microeco object from CPPM Combined.RDS


# Filter to 10% prevalence (Taxa present in 10% of samples)
ps <- meco2phyloseq(m16)

library(metagMisc)
psf <- phyloseq_filter_prevalence(ps, prev.trh = 0.1, abund.trh = NULL, threshold_condition = "OR", abund.type = "total")
psf

m16 <- phyloseq2meco(psf)

# Compare network modules 
# Fusarium with brome and fusarium without brome


# by location
bz_wfb <- clone(m16)
bz_nfb <- clone(m16)
sarc_wfb <- clone(m16)
sarc_nfb <- clone(m16)

bz_wfb$sample_table %<>% subset(Location == "Bozeman" & B.tectorum == "B. tectorum" & F..pseudograminearumPresence == "F. pseudograminearum")
bz_nfb$sample_table %<>% subset(Location == "Bozeman" & B.tectorum == "No B. tectorume" & F..pseudograminearumPresence == "No F. pseudograminearum")
sarc_wfb$sample_table %<>% subset(Location == "Huntley" & B.tectorum == "B. tectorum" & F..pseudograminearumPresence == "F. pseudograminearum")
sarc_nfb$sample_table %<>% subset(Location == "Huntley" & B.tectorum == "No B. tectorum" & F..pseudograminearumPresence == "No F. pseudograminearum")

bz_wfb$tidy_dataset()
bz_nfb$tidy_dataset()
sarc_wfb$tidy_dataset()
sarc_nfb$tidy_dataset()

# Construct Networks
bz_wfb_net <- trans_network$new(dataset = bz_wfb, cor_method = NULL, taxa_level = "OTU", filter_thres = 0.0002)
bz_wfb_net$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb", icov.select.params(ncores = 40), 
                       add_taxa_name = c("Kingdom", "Phylum", "Family", "Genus"), verbose=TRUE)
saveRDS(bz_wfb_net, "bz_wfb_net.RDS")

bz_nfb_net <- trans_network$new(dataset = bz_nfb, cor_method = NULL, taxa_level = "OTU", filter_thres = 0.0002)
bz_nfb_net$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb", icov.select.params(ncores = 40), 
                       add_taxa_name = c("Kingdom", "Phylum", "Family", "Genus"), verbose=TRUE)
saveRDS(bz_nfb_net, "bz_nfb_net.RDS")

sarc_wfb_net <- trans_network$new(dataset = sarc_wfb, cor_method = NULL, taxa_level = "OTU", filter_thres = 0.0002)
sarc_wfb_net$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb", icov.select.params(ncores = 40), 
                       add_taxa_name = c("Kingdom", "Phylum", "Family", "Genus"), verbose=TRUE)
saveRDS(sarc_wfb_net, "sarc_wfb_net.RDS")

sarc_nfb_net <- trans_network$new(dataset = sarc_nfb, cor_method = NULL, taxa_level = "OTU", filter_thres = 0.0002)
sarc_nfb_net$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb", icov.select.params(ncores = 40), 
                       add_taxa_name = c("Kingdom", "Phylum", "Family", "Genus"), verbose=TRUE)
saveRDS(sarc_nfb_net, "sarc_nfb_net.RDS")

# Import previously constructed networks
bz_wfb_net <- readRDS("bz_wfb_net.RDS")
bz_nfb_net <- readRDS("bz_nfb_net.RDS")
sarc_wfb_net <- readRDS("sarc_wfb_net.RDS")
sarc_nfb_net <- readRDS("sarc_nfb_net.RDS")

bz_wfb_net$cal_module(method = "cluster_fast_greedy")
bz_nfb_net$cal_module(method = "cluster_fast_greedy")
sarc_wfb_net$cal_module(method = "cluster_fast_greedy")
sarc_nfb_net$cal_module(method = "cluster_fast_greedy")

bz_wfb_net$save_network(filepath = "bz_wfb_net.gexf")
bz_nfb_net$save_network(filepath = "bz_nfb_net.gexf")
sarc_wfb_net$save_network(filepath = "sarc_wfb_net.gexf")
sarc_nfb_net$save_network(filepath = "sarc_nfb_net.gexf")

# View directly in Cytoscape
# https://cytoscape.org/cytoscape-automation/for-scripters/R/notebooks/Network-functions-and-visualization.nb.html
library(RCy3)

bz_wfb_igraph_net <- bz_wfb_net$res_network
bz_nfb_igraph_net <- bz_nfb_net$res_network
sarc_wfb_igraph_net <- sarc_wfb_net$res_network
sarc_nfb_igraph_net <- sarc_nfb_net$res_network

bz_wfb_degAll <- igraph::degree(bz_wfb_igraph_net, v = igraph::V(bz_wfb_igraph_net), mode = "all")
bz_wfb_betAll <- igraph::betweenness(bz_wfb_igraph_net, v = igraph::V(bz_wfb_igraph_net), directed = FALSE) / (((igraph::vcount(bz_wfb_igraph_net) - 1) * (igraph::vcount(bz_wfb_igraph_net)-2)) / 2)
bz_wfb_betAll.norm <- (bz_wfb_betAll - min(bz_wfb_betAll))/(max(bz_wfb_betAll) - min(bz_wfb_betAll))
bz_wfb_dsAll <- igraph::similarity(method = "dice", bz_wfb_igraph_net, vids = igraph::V(bz_wfb_igraph_net), mode = "all")
bz_wfb_igraph_net <- igraph::set_vertex_attr(bz_wfb_igraph_net, "degree", index = igraph::V(bz_wfb_igraph_net), value = bz_wfb_degAll)
bz_wfb_igraph_net <- igraph::set_vertex_attr(bz_wfb_igraph_net, "betweenness", index = igraph::V(bz_wfb_igraph_net), value = bz_wfb_betAll.norm)
bz_wfb_igraph_net <- igraph::set_edge_attr(bz_wfb_igraph_net, "weight", index = igraph::E(bz_wfb_igraph_net), value = NULL)

bz_nfb_degAll <- igraph::degree(bz_nfb_igraph_net, v = igraph::V(bz_nfb_igraph_net), mode = "all")
bz_nfb_betAll <- igraph::betweenness(bz_nfb_igraph_net, v = igraph::V(bz_nfb_igraph_net), directed = FALSE) / (((igraph::vcount(bz_nfb_igraph_net) - 1) * (igraph::vcount(bz_nfb_igraph_net)-2)) / 2)
bz_nfb_betAll.norm <- (bz_nfb_betAll - min(bz_nfb_betAll))/(max(bz_nfb_betAll) - min(bz_nfb_betAll))
bz_nfb_dsAll <- igraph::similarity(method = "dice", bz_nfb_igraph_net, vids = igraph::V(bz_nfb_igraph_net), mode = "all")
bz_nfb_igraph_net <- igraph::set_vertex_attr(bz_nfb_igraph_net, "degree", index = igraph::V(bz_nfb_igraph_net), value = bz_nfb_degAll)
bz_nfb_igraph_net <- igraph::set_vertex_attr(bz_nfb_igraph_net, "betweenness", index = igraph::V(bz_nfb_igraph_net), value = bz_nfb_betAll.norm)
bz_nfb_igraph_net <- igraph::set_edge_attr(bz_nfb_igraph_net, "weight", index = igraph::E(bz_nfb_igraph_net), value = NULL)

sarc_wfb_degAll <- igraph::degree(sarc_wfb_igraph_net, v = igraph::V(sarc_wfb_igraph_net), mode = "all")
sarc_wfb_betAll <- igraph::betweenness(sarc_wfb_igraph_net, v = igraph::V(sarc_wfb_igraph_net), directed = FALSE) / (((igraph::vcount(sarc_wfb_igraph_net) - 1) * (igraph::vcount(sarc_wfb_igraph_net)-2)) / 2)
sarc_wfb_betAll.norm <- (sarc_wfb_betAll - min(sarc_wfb_betAll))/(max(sarc_wfb_betAll) - min(sarc_wfb_betAll))
sarc_wfb_dsAll <- igraph::similarity(method = "dice", sarc_wfb_igraph_net, vids = igraph::V(sarc_wfb_igraph_net), mode = "all")
sarc_wfb_igraph_net <- igraph::set_vertex_attr(sarc_wfb_igraph_net, "degree", index = igraph::V(sarc_wfb_igraph_net), value = sarc_wfb_degAll)
sarc_wfb_igraph_net <- igraph::set_vertex_attr(sarc_wfb_igraph_net, "betweenness", index = igraph::V(sarc_wfb_igraph_net), value = sarc_wfb_betAll.norm)
sarc_wfb_igraph_net <- igraph::set_edge_attr(sarc_wfb_igraph_net, "weight", index = igraph::E(sarc_wfb_igraph_net), value = NULL)

sarc_nfb_degAll <- igraph::degree(sarc_nfb_igraph_net, v = igraph::V(sarc_nfb_igraph_net), mode = "all")
sarc_nfb_betAll <- igraph::betweenness(sarc_nfb_igraph_net, v = igraph::V(sarc_nfb_igraph_net), directed = FALSE) / (((igraph::vcount(sarc_nfb_igraph_net) - 1) * (igraph::vcount(sarc_nfb_igraph_net)-2)) / 2)
sarc_nfb_betAll.norm <- (sarc_nfb_betAll - min(sarc_nfb_betAll))/(max(sarc_nfb_betAll) - min(sarc_nfb_betAll))
sarc_nfb_dsAll <- igraph::similarity(method = "dice", sarc_nfb_igraph_net, vids = igraph::V(sarc_nfb_igraph_net), mode = "all")
sarc_nfb_igraph_net <- igraph::set_vertex_attr(sarc_nfb_igraph_net, "degree", index = igraph::V(sarc_nfb_igraph_net), value = sarc_nfb_degAll)
sarc_nfb_igraph_net <- igraph::set_vertex_attr(sarc_nfb_igraph_net, "betweenness", index = igraph::V(sarc_nfb_igraph_net), value = sarc_nfb_betAll.norm)
sarc_nfb_igraph_net <- igraph::set_edge_attr(sarc_nfb_igraph_net, "weight", index = igraph::E(sarc_nfb_igraph_net), value = NULL)


# Create network in Cytoscape
createNetworkFromIgraph(bz_wfb_igraph_net, new.title='Bozeman Fusarium Brome')
createNetworkFromIgraph(bz_nfb_igraph_net, new.title='Bozeman No Fusarium Brome')
createNetworkFromIgraph(sarc_wfb_igraph_net, new.title='Huntley Fusarium Brome')
createNetworkFromIgraph(sarc_nfb_igraph_net, new.title='Huntley No Fusarium Brome')

# Import network tables from Cytoscape
bz_wfb_at <- read.csv("Bozeman Fusarium Brome default node.csv")
bz_nfb_at <- read.csv("Bozeman No Fusarium Brome default node.csv")
sarc_wfb_at <- read.csv("Huntley Fusarium Brome default node.csv")
sarc_nfb_at <- read.csv("Huntley No Fusarium Brome default node.csv")

bz_wfb_at$Treatment <- "bz_wfb"
bz_nfb_at$Treatment <- "bz_nfb"
sarc_wfb_at$Treatment <- "sarc_wfb"
sarc_nfb_at$Treatment <- "sarc_nfb"

mdf <- rbind(bz_wfb_at, bz_nfb_at, sarc_wfb_at, sarc_nfb_at)

# Compare to Gephi
bz_wfb_at1 <- read.csv("bz_wfb_table_gephi.csv")
bz_nfb_at1 <- read.csv("bz_nfb_table_gephi.csv")
sarc_wfb_at1 <- read.csv("sarc_wfb_table_gephi.csv")
sarc_nfb_at1 <- read.csv("bz_nfb_table_gephi.csv")

bz_wfb_at1$Treatment <- "bz_wfb"
bz_nfb_at1$Treatment <- "bz_nfb"
sarc_wfb_at1$Treatment <- "sarc_wfb"
sarc_nfb_at1$Treatment <- "sarc_nfb"

colnames(sarc_wfb_at1)
bz_wfb_at1 <- bz_wfb_at1[,-c(19,22)]
sarc_wfb_at1 <- sarc_wfb_at1[,-c(19,22)]

mdf1 <- rbind(bz_wfb_at1, bz_nfb_at1, sarc_wfb_at1, sarc_nfb_at1)

# Plots
library(tidyr)
library(RColorBrewer)
library(ggpubr)
library(phylosmith)

df_melt <- mdf[,c(1,3:5,26,27)] # Subset to data to compare.

df_long <- df_melt %>% 
  pivot_longer(cols = 1:5, names_to = "Variable", values_to = "Value")

my_comp <- list(c("bz_wfb", "bz_nfb"), c("sarc_wfb", "sarc_nfb"))
v1 <- ggviolin(df_long, x = "Treatment", y = "Value", fill = "Treatment", orientation = "horizontal", palette = brewer.pal(4,"Dark2"),
               add = "boxplot", add.params = list(fill = "white")) + 
  stat_compare_means(comparisons = my_comp, label = "p.signif")+ # Add significance levels
  facet_wrap(~Variable, scales = "free") + theme(legend.position = "none") +
  theme(plot.margin = margin(t = 0, r = 0.5, b = 0, l = 0, unit = "cm")) 
v1 

df_melt1 <- mdf1[,c(11,13:15,22)] # Subset to data to compare.

df_long1 <- df_melt1 %>% 
  pivot_longer(cols = 1:4, names_to = "Variable", values_to = "Value")

my_comp1 <- list(c("bz_wfb", "bz_nfb"), c("sarc_wfb", "sarc_nfb"))
v2 <- ggviolin(df_long1, x = "Treatment", y = "Value", fill = "Treatment", orientation = "horizontal", palette = brewer.pal(4,"Dark2"),
               add = "boxplot", add.params = list(fill = "white")) + 
  stat_compare_means(comparisons = my_comp, label = "p.signif")+ # Add significance levels
  facet_wrap(~Variable, scales = "free") + theme(legend.position = "none") +
  theme(plot.margin = margin(t = 0, r = 0.5, b = 0, l = 0, unit = "cm")) 
v2 

# Hub taxa calculated with Cytohubba in Cytoscape and .csvs imported into R for analysis
bzb_hub <- read.csv("Bozeman Fusarium Brome_MCC_top10 default node.csv")
bznb_hub <- read.csv("Bozeman No Fusarium Brome_MCC_top10 default node.csv")
sarcb_hub <- read.csv("Huntley Fusarium Brome_MCC_top10 default node.csv")
sarcnb_hub <- read.csv("Huntley No Fusarium Brome_MCC_top10 default node.csv")

bzb_hub$Location <- "Bozeman"
bznb_hub$Location <- "Bozeman"
sarcb_hub$Location <- "Huntley"
sarcnb_hub$Location <- "Huntley"

bzb_hub$Fusarium <- "Fusarium"
bznb_hub$Fusarium <- "No Fusarium"
sarcb_hub$Fusarium <- "Fusarium"
sarcnb_hub$Fusarium <- "No Fusarium"

hub_df <- rbind(bzb_hub, bznb_hub, sarcb_hub, sarcnb_hub)
write.csv(hub_df, "Hub Taxa Loc Fusarium.csv")

unique(hub_df$Phylum)
hub_df$Phylum <- factor(hub_df$Phylum, levels = c("Acidobacteria", "Actinobacteria",  "Bacteroidetes", "Firmicutes", "Armatimonadetes", "Proteobacteria", "Cyanobacteria"))
ggbarplot(hub_df, x = "Fusarium", y= "RelativeAbundance", fill = "Phylum", facet.by = "Location") + ylab("Relative Abundance (%)")

unique(hub_df$Family)
ggbarplot(hub_df, x = "Fusarium", y= "RelativeAbundance", fill = "Family", facet.by = "Location") + ylab("Relative Abundance (%)")

# Venn Diagram
library(ggvenn)

# Phylum
x2 <- list('Bozeman Fusarium' = bzb_hub$Phylum, 'Bozeman No Fusarium' = bznb_hub$Phylum, "Huntley No Fusarium" = sarcnb_hub$Phylum, 'Huntley Fusarium' = sarcb_hub$Phylum)

pdf("Figures/Venn Phylum Hub.pdf", width = 10, height = 8)
ggvenn(
  x2, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#CD80C2FF" ,"#ED73C2FF"),
  stroke_size = 0, set_name_size = 4, show_outside = "always"
)
dev.off() 

# Family
x3 <- list('Bozeman Fusarium' = bzb_hub$Family, 'Bozeman No Fusarium' = bznb_hub$Family, "Huntley No Fusarium" = sarcnb_hub$Family, 'Huntley Fusarium' = sarcb_hub$Family)

pdf("Figures/Venn Family Hub.pdf", width = 10, height = 8)
ggvenn(
  x3, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#CD80C2FF" ,"#ED73C2FF"),
  stroke_size = 0, set_name_size = 4, show_outside = "always"
)
dev.off() 

# Genus
x4 <- list('Bozeman Fusarium' = bzb_hub$Genus, 'Bozeman No Fusarium' = bznb_hub$Genus, "Huntley No Fusarium" = sarcnb_hub$Genus, 'Huntley Fusarium' = sarcb_hub$Genus)

pdf("Figures/Venn Genus Hub.pdf", width = 10, height = 8)
ggvenn(
  x4, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#CD80C2FF" ,"#ED73C2FF"),
  stroke_size = 0, set_name_size = 4, show_outside = "always"
)
dev.off() 


########## END By Location ##############


##################
# Fusarium Brome #
##################

# Across location years comparing Fusarium and Brome to No Fusarium or Brome
wfb <- clone(m16)
nfb <- clone(m16)
wfnb <- clone(m16)
nfwb <- clone(m16)

wfb$sample_table %<>% subset(Weeds == "Brome" & FusariumPresence == "Fusarium")
nfb$sample_table %<>% subset(Weeds == "No Brome" & FusariumPresence == "No Fusarium")
wfnb$sample_table %<>% subset(Weeds == "No Brome" & FusariumPresence == "Fusarium")
nfwb$sample_table %<>% subset(Weeds == "Brome" & FusariumPresence == "No Fusarium")

wfb$tidy_dataset()
nfb$tidy_dataset()
wfnb$tidy_dataset()
nfwb$tidy_dataset()

# Construct Networks
wfb_net <- trans_network$new(dataset = wfb, cor_method = NULL, taxa_level = "OTU", filter_thres = 0.0005)
wfb_net$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb", icov.select.params(ncores = 40), 
                       add_taxa_name = c("Kingdom", "Phylum", "Family", "Genus"), verbose=TRUE)
saveRDS(wfb_net, "wfb_net.RDS")

nfb_net <- trans_network$new(dataset = nfb, cor_method = NULL, taxa_level = "OTU", filter_thres = 0.0005)
nfb_net$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb", icov.select.params(ncores = 40), 
                       add_taxa_name = c("Kingdom", "Phylum", "Family", "Genus"), verbose=TRUE)
saveRDS(nfb_net, "nfb_net.RDS")

wfnb_net <- trans_network$new(dataset = wfnb, cor_method = NULL, taxa_level = "OTU", filter_thres = 0.0005)
wfnb_net$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb", icov.select.params(ncores = 40), 
                         add_taxa_name = c("Kingdom", "Phylum", "Family", "Genus"), verbose=TRUE)
saveRDS(wfnb_net, "wfnb_net.RDS")

nfwb_net <- trans_network$new(dataset = nfwb, cor_method = NULL, taxa_level = "OTU", filter_thres = 0.0005)
nfwb_net$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb", icov.select.params(ncores = 40), 
                         add_taxa_name = c("Kingdom", "Phylum", "Family", "Genus"), verbose=TRUE)
saveRDS(nfwb_net, "nfwb_net.RDS")


# Import previously constructed networks
wfb_net <- readRDS("wfb_net.RDS")
nfb_net <- readRDS("nfb_net.RDS")
wfnb_net <- readRDS("wfnb_net.RDS")
nfwb_net <- readRDS("nfwb_net.RDS")

wfb_net$cal_module(method = "cluster_fast_greedy")
nfb_net$cal_module(method = "cluster_fast_greedy")
wfnb_net$cal_module(method = "cluster_fast_greedy")
nfwb_net$cal_module(method = "cluster_fast_greedy")

wfb_net$save_network(filepath = "wfb_net.gexf")
nfb_net$save_network(filepath = "nfb_net.gexf")
wfnb_net$save_network(filepath = "wfnb_net.gexf")
nfwb_net$save_network(filepath = "nfwb_net.gexf")

wfb_net$plot_taxa_roles(use_type = 1)
nfb_net$plot_taxa_roles(use_type = 1)
wfnb_net$plot_taxa_roles(use_type = 1)
nfwb_net$plot_taxa_roles(use_type = 1)

wfb_net$cal_eigen()
nfb_net$cal_eigen()
wfnb_net$cal_eigen()
nfwb_net$cal_eigen()
wfb$sample_table
wfb_env <- trans_env$new(dataset = wfb, add_data = wfb$sample_table[,c(12:16)])
nfb_env <- trans_env$new(dataset = nfb, add_data = nfb$sample_table[,c(12:16)])
wfnb_env <- trans_env$new(dataset = wfnb, add_data = wfnb$sample_table[,c(12:16)])
nfwb_env <- trans_env$new(dataset = nfwb, add_data = nfwb$sample_table[,c(12:16)])

# calculate correlations
wfb_env$cal_cor(add_abund_table = wfb_net$res_eigen)
nfb_env$cal_cor(add_abund_table = nfb_net$res_eigen)
wfnb_env$cal_cor(add_abund_table = wfnb_net$res_eigen)
nfwb_env$cal_cor(add_abund_table = nfwb_net$res_eigen)

# plot the correlation heatmap
wfb_env$plot_cor()
nfb_env$plot_cor()
wfnb_env$plot_cor()
nfwb_env$plot_cor()

# Compare networks
library(meconetcomp)

net_list <- list()
net_list$nfb <- nfb_net
net_list$wfb <- wfb_net
net_list$wfnb <- wfnb_net
net_list$nfwb <- nfwb_net

net_list %<>% cal_module(undirected_method = "cluster_fast_greedy")
# nfb: 17 modules, wfb: 13 modules, wfnb: 15 modules, nfwb: 18 modules
net_list %<>% get_node_table(node_roles = TRUE) %>% get_edge_table

# Subset Networks
sn <- edge_comp(net_list)
sub_net <- trans_venn$new(sn)
# convert intersection result to a microtable object
sub_net1 <- sub_net$trans_comm()
# extract the intersection of all four networks 
# please use colnames(tmp2$otu_table) to find the required name
Intersec_all <- subset_network(net_list, venn = sub_net1, name = "nfb&wfb&wfnb&nfwb")
# Intersec_all is a trans_network object
# for example, save Intersec_all as gexf format

Intersec_all$save_network("Intersect_net.gexf")

# Module Comparison
wfb_at <- read.csv("Fusarium Brome default node.csv")
nfb_at <- read.csv("No Fusarium Brome default node.csv")
wfnb_at <- read.csv("Fusarium No Brome default node.csv")
nfwb_at <- read.csv("No Fusarium Brome default node.csv")

nfb_M1 <- subset(wfb_at, module == "M1")
wfb_M1 <- subset(nfb_at, module == "M1")
wfnb_M1 <- subset(wfnb_at, module == "M1")
nfwb_M1 <- subset(nfwb_at, module == "M1")

M1 <- rbind(nfb_M1, wfb_M1, wfnb_M1, nfwb_M1)

M1_melt <- M1[,c(1:3,5,19,26)] # Subset to data to compare.

M1_long1 <- M1_melt %>% 
  pivot_longer(cols = 1:6, names_to = "Variable", values_to = "Value")

m1_comp <- list(c("nfb", "wfb"), c("nfb", "nfwb"), c("nfb", "wfnb"))

M1p <- ggviolin(M1_long1, x = "Treatment", y = "Value", fill = "Treatment", orientation = "horizontal", palette = brewer.pal(4,"Dark2"),
               add = "boxplot", add.params = list(fill = "white")) + 
  stat_compare_means(comparisons = m1_comp, label = "p.signif")+ # Add significance levels
  facet_wrap(~Variable, scales = "free") + theme(legend.position = "none") +
  theme(plot.margin = margin(t = 0, r = 0.5, b = 0, l = 0, unit = "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
M1p 

# Module M2
nfb_M2 <- subset(wfb_at, module == "M2")
wfb_M2 <- subset(nfb_at, module == "M2")
wfnb_M2 <- subset(wfnb_at, module == "M2")
nfwb_M2 <- subset(nfwb_at, module == "M2")

M2 <- rbind(nfb_M2, wfb_M2, wfnb_M2, nfwb_M2)

M2_melt <- M2[,c(1:3,5,19,26)] # Subset to data to compare.

M2_long1 <- M2_melt %>% 
  pivot_longer(cols = 1:6, names_to = "Variable", values_to = "Value")

M2p <- ggviolin(M2_long1, x = "Treatment", y = "Value", fill = "Treatment", orientation = "horizontal", palette = brewer.pal(4,"Dark2"),
                add = "boxplot", add.params = list(fill = "white")) + 
  stat_compare_means(comparisons = m1_comp, label = "p.signif")+ # Add significance levels
  facet_wrap(~Variable, scales = "free") + theme(legend.position = "none") +
  theme(plot.margin = margin(t = 0, r = 0.5, b = 0, l = 0, unit = "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
M2p 

# Module M3
nfb_M3 <- subset(wfb_at, module == "M3")
wfb_M3 <- subset(nfb_at, module == "M3")
wfnb_M3 <- subset(wfnb_at, module == "M3")
nfwb_M3 <- subset(nfwb_at, module == "M3")

M3 <- rbind(nfb_M3, wfb_M3, wfnb_M3, nfwb_M3)

M3_melt <- M3[,c(1:3,5,19,26,27)] # Subset to data to compare.

M3_long1 <- M3_melt %>% 
  pivot_longer(cols = 1:6, names_to = "Variable", values_to = "Value")

M3p <- ggviolin(M3_long1, x = "Treatment", y = "Value", fill = "Treatment", orientation = "horizontal", palette = brewer.pal(4,"Dark2"),
                add = "boxplot", add.params = list(fill = "white")) + 
  stat_compare_means(comparisons = m1_comp, label = "p.signif")+ # Add significance levels
  facet_wrap(~Variable, scales = "free") + theme(legend.position = "none") +
  theme(plot.margin = margin(t = 0, r = 0.5, b = 0, l = 0, unit = "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
M3p 

# Module M4
nfb_M4 <- subset(wfb_at, module == "M4")
wfb_M4 <- subset(nfb_at, module == "M4")
wfnb_M4 <- subset(wfnb_at, module == "M4")
nfwb_M4 <- subset(nfwb_at, module == "M4")

M4 <- rbind(nfb_M4, wfb_M4, wfnb_M4, nfwb_M4)

M4_melt <- M4[,c(1:3,5,19,26,27)] # Subset to data to compare.

M4_long1 <- M4_melt %>% 
  pivot_longer(cols = 1:6, names_to = "Variable", values_to = "Value")

M4p <- ggviolin(M4_long1, x = "Treatment", y = "Value", fill = "Treatment", orientation = "horizontal", palette = brewer.pal(4,"Dark2"),
                add = "boxplot", add.params = list(fill = "white")) + 
  stat_compare_means(comparisons = m1_comp, label = "p.signif")+ # Add significance levels
  facet_wrap(~Variable, scales = "free") + theme(legend.position = "none") +
  theme(plot.margin = margin(t = 0, r = 0.5, b = 0, l = 0, unit = "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
M4p 

ggarrange(M1p, M2p, M3p, M4p, nrow = 2, ncol = 2, labels = c("a", "b", "c", "d"))


module_list <- list()
module_list$`No Fusarium \n or Brome` <- nfb_module$Phylum
module_list$`With Fusarium \n and Brome` <- wfb_module$Phylum
module_list$`With Fusarium \n No Brome` <- wfnb_module$Phylum
module_list$`No Fusarium \n With Brome` <- nfwb_module$Phylum

library(ggvenn)
mv1 <-ggvenn(module_list, fill_color = c("#1B9E77", "#D95F02", "#7570B3", "#66A61E"), stroke_color = NA, label_sep = "\n", padding = 0.2)
mv1
# 45.5% of phyla in module 1 are shared by all treatments

pdf("Figures/M1 Shared Phyla.pdf", width = 8, height = 8)
ggvenn(module_list, fill_color = c("#1B9E77", "#D95F02", "#7570B3", "#66A61E"), stroke_color = NA, padding = 0.2) 
dev.off()

wfb_module$Treatment <- "With Fusarium and Brome"
nfb_module$Treatment <- "No Fusarium No Brome"
wfnb_module$Treatment <- "With Fusarium No Brome"
nfwb_module$Treatment <- "No Fusarium With Brome"

M1 <- rbind(wfb_module, nfb_module, wfnb_module, nfwb_module)
M1 <- M1[,-c(1:10)]
write.csv(M1, "Module M1 Taxa.csv")

# Compare phylogenetic distances of paired nodes in edges
# filter useless features to speed up the calculation

node_names <- unique(unlist(lapply(net_list, function(x){colnames(x$data_abund)})))
filter_amp <- microeco::clone(m16)
filter_amp$otu_table <- filter_amp$otu_table[node_names, ]
filter_amp$tidy_dataset()
# obtain phylogenetic distance matrix
phylogenetic_distance <- as.matrix(cophenetic(filter_amp$phylo_tree))
# use both the positive and negative labels
tmp <- edge_node_distance$new(network_list = net_list, dis_matrix = phylogenetic_distance, label = c("+", "-"))
tmp$cal_diff(method = "anova")

# show different modules with at least 10 nodes and positive edges
tmp <- edge_node_distance$new(network_list = net_list, dis_matrix = phylogenetic_distance, 
                              label = "+", with_module = TRUE, module_thres = 25)
tmp$cal_diff(method = "anova")
g1 <- tmp$plot(add = "none", add_sig = TRUE, add_sig_text_size = 5) + ylab("Phylogenetic distance")
g1
ggsave("Figures/amp_phylo_distance_modules.pdf", g1, width = 8, height = 6)

amp_network_edgetax <- edge_tax_comp(net_list, taxrank = "Phylum", label = "+", rel = TRUE)
# filter the features with small number
amp_network_edgetax <- amp_network_edgetax[apply(amp_network_edgetax, 1, mean) > 0.01, ]
# visualization
g2 <- pheatmap::pheatmap(amp_network_edgetax, display_numbers = TRUE)
g2

ggsave("soil_amp_edge_tax_comp.pdf", g1, width = 7, height = 7)

# calculate global properties of all sub-networks
tmp <- subnet_property(net_list)
# then prepare the data for the correlation analysis
# use sample names (second column) as rownames
rownames(tmp) <- tmp[, 2]
# delete first two columns (network name and sample name)
tmp <- tmp[, -c(1:2)]
# load ready-made abiotic factor and diversity table
tmp1 <- trans_env$new(dataset = m16, add_data = m16$sample_table[, c(12:16,19)])
tmp1$cal_cor(use_data = "other", by_group = "Fus_Brome", add_abund_table = tmp, method = "spearman")
# generate correlation heatmap
g3 <- tmp1$plot_cor()
g3
ggsave("Figures/soil_amp_subnet_property.pdf", g1, width = 11, height = 5)

# Robustness of Networks
rb <- robustness$new(net_list, remove_strategy = c("edge_rand", "edge_strong", "node_rand", "node_degree_high"), 
                      remove_ratio = seq(0, 0.99, 0.1), measure = c("Eff", "Eigen", "Pcr"), run = 10)
#View(rb$res_table)
#View(rb$res_summary)
rb$plot(linewidth = 1)
ggsave("Figures/Network Robustness.pdf", rb$plot(linewidth = 1), width = 10, height = 6)

ex1 <- rb$res_table %>% .[.$remove_strategy == "node_rand" & .$measure == "Eigen", ]

t1 <- trans_env$new(dataset = NULL, add_data = ex1)
t1$dataset$sample_table <- t1$data_env
t1$plot_scatterfit(x = "remove_ratio", y = "value", type = "cor", group = "Network") + 
  xlab("Ratio of randomly removed nodes") + ylab("Network connectivity") + theme(axis.title = element_text(size = 15))

vul_table <- vulnerability(net_list)
View(vul_table)
ggboxplot(vul_table, x= "Network", y="vulnerability")

c1 <- cohesionclass$new(net_list)
View(c1$res_list$sample)
View(c1$res_list$feature)
c1$cal_diff(method = "anova", measure = "c_pos")

cp1 <- c1$plot(measure = "c_pos") + ylab("Cohesion") + font("ylab", size = 12) +
  theme(plot.margin = margin(t = 0, r = 0.5, b = 0, l = 0, unit = "cm")) 
cp1

pdf("Figures/Cohesion.pdf", width = 10, height = 6)
c1$plot(measure = "r_pos") + ylab("Cohesion")
dev.off()

r1 <- cohesionclass$new(net_list)
r1$cal_diff(method = "anova", measure = "r_pos")
r1$plot("r_pos")

nc1 <- cohesionclass$new(net_list)
nc1$cal_diff(method = "anova", measure = "c_neg")
nc1$plot("c_neg")

nr1 <- cohesionclass$new(net_list)
nr1$cal_diff(method = "anova", measure = "r_neg")
nr1$plot("r_neg")

cdf <- c1$res_list$sample
rdf <- c1$res_list$feature

cdf$Total <- abs(cdf$c_neg)+cdf$c_pos
rdf$Total <- abs(rdf$r_neg)+rdf$r_pos

ggboxplot(cdf, x = "network", y="Total") + ylim(0,0.5)

rdf$network <- factor(rdf$network, levels = c("nfb", "nfwb", "wfb", "wfnb"))
aov_result <- aov(Total ~ network, data = rdf)
tukey_result <- TukeyHSD(aov_result)

cld_result <- multcompLetters4(aov_result, tukey_result)

# Extract letters and add to data
letters_df <- data.frame(
  group = names(cld_result$network$Letters),
  letter = cld_result$network$Letters
)


cp1 <- ggviolin(rdf, x = "network", y="Total", fill = "network", palette = "Dark2", add = "boxplot", width = 0.5, add.params = list(fill = "white")) + ylab("Total Cohesion") + 
  xlab("Treatment") + theme(legend.position = "none", axis.text.x = element_blank()) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  geom_text(data = letters_df, aes(x = group, y = Inf, label = letter),
            vjust = 1.5, size = 4, inherit.aes = FALSE) + theme(plot.margin = margin(t = 0, r = 0.5, b = 0, l = 0, unit = "cm")) 
cp1

mean(rdf$network=="nfb") # 0.25
mean(rdf$network=="nfwb") # 0.25
mean(rdf$network=="wfb")
mean(rdf$network=="wfnb")

# Compare nodes across networks
# obtain the node distributions by searching the res_node_table in the object
n_comp <- node_comp(net_list, property = "name")
tmp <- cal_network_attr(net_list)
tmp

net_list %<>% cal_module(undirected_method = "cluster_fast_greedy")
net_list %<>% get_node_table(node_roles = TRUE) %>% get_edge_table

# Rename
n_comp$sample_table

rownames(n_comp$sample_table)[rownames(n_comp$sample_table) == "nfb"] <- "No F. pseudograminearum No B. tectorum"
rownames(n_comp$sample_table)[rownames(n_comp$sample_table) == "wfb"] <- "F. pseudograminearum B. tectorum"
rownames(n_comp$sample_table)[rownames(n_comp$sample_table) == "wfnb"] <- "F. pseudograminearum No B. tectorum"
rownames(n_comp$sample_table)[rownames(n_comp$sample_table) == "nfwb"] <- "No F. pseudograminearum B. tectorum"

n_comp$otu_table <- n_comp$otu_table %>% dplyr::rename("No F. pseudograminearum No B. tectorum" = "nfb", "F. pseudograminearum B. tectorum" = "wfb",
                                                "F. pseudograminearum No B. tectorum" = "wfnb", "No F. pseudograminearum B. tectorum" = "nfwb")

n_comp$otu_table <- n_comp$otu_table[,c(2,3,4,1)]
n_comp$sample_table <- n_comp$sample_table[c(2,3,4,1),]

# obtain nodes intersection
vcomp <- trans_venn$new(n_comp, ratio = "numratio")
vcomp$plot_venn(fill_color = TRUE)
g1 <- vcomp$plot_venn(fill_color = FALSE)
g1
ggsave("Figures/FB_node_overlap.pdf", g1, width = 7, height = 6)

np <- vcomp$plot_bar(left_plot = TRUE, bottom_height = 0.5, left_width = 0.25, up_bar_fill = "grey50", left_bar_fill = "grey50", bottom_point_color = "black", left_x_text_size = 9)
np


# calculate jaccard distance to reflect the overall differences of networks
n_comp$cal_betadiv(method = "bray")
n_comp$beta_diversity$bray


# View directly in Cytoscape
# https://cytoscape.org/cytoscape-automation/for-scripters/R/notebooks/Network-functions-and-visualization.nb.html
library(RCy3)

wfb_igraph_net <- wfb_net$res_network
nfb_igraph_net <- nfb_net$res_network
wfnb_igraph_net <- wfnb_net$res_network
nfwb_igraph_net <- nfwb_net$res_network

wfb_degAll <- igraph::degree(wfb_igraph_net, v = igraph::V(wfb_igraph_net), mode = "all")
wfb_betAll <- igraph::betweenness(wfb_igraph_net, v = igraph::V(wfb_igraph_net), directed = FALSE) / (((igraph::vcount(wfb_igraph_net) - 1) * (igraph::vcount(wfb_igraph_net)-2)) / 2)
wfb_betAll.norm <- (wfb_betAll - min(wfb_betAll))/(max(wfb_betAll) - min(wfb_betAll))
wfb_dsAll <- igraph::similarity(method = "dice", wfb_igraph_net, vids = igraph::V(wfb_igraph_net), mode = "all")
wfb_igraph_net <- igraph::set_vertex_attr(wfb_igraph_net, "degree", index = igraph::V(wfb_igraph_net), value = wfb_degAll)
wfb_igraph_net <- igraph::set_vertex_attr(wfb_igraph_net, "betweenness", index = igraph::V(wfb_igraph_net), value = wfb_betAll.norm)
wfb_igraph_net <- igraph::set_edge_attr(wfb_igraph_net, "weight", index = igraph::E(wfb_igraph_net), value = NULL)

nfb_degAll <- igraph::degree(nfb_igraph_net, v = igraph::V(nfb_igraph_net), mode = "all")
nfb_betAll <- igraph::betweenness(nfb_igraph_net, v = igraph::V(nfb_igraph_net), directed = FALSE) / (((igraph::vcount(nfb_igraph_net) - 1) * (igraph::vcount(nfb_igraph_net)-2)) / 2)
nfb_betAll.norm <- (nfb_betAll - min(nfb_betAll))/(max(nfb_betAll) - min(nfb_betAll))
nfb_dsAll <- igraph::similarity(method = "dice", nfb_igraph_net, vids = igraph::V(nfb_igraph_net), mode = "all")
nfb_igraph_net <- igraph::set_vertex_attr(nfb_igraph_net, "degree", index = igraph::V(nfb_igraph_net), value = nfb_degAll)
nfb_igraph_net <- igraph::set_vertex_attr(nfb_igraph_net, "betweenness", index = igraph::V(nfb_igraph_net), value = nfb_betAll.norm)
nfb_igraph_net <- igraph::set_edge_attr(nfb_igraph_net, "weight", index = igraph::E(nfb_igraph_net), value = NULL)

wfnb_degAll <- igraph::degree(wfnb_igraph_net, v = igraph::V(wfnb_igraph_net), mode = "all")
wfnb_betAll <- igraph::betweenness(wfnb_igraph_net, v = igraph::V(wfnb_igraph_net), directed = FALSE) / (((igraph::vcount(wfnb_igraph_net) - 1) * (igraph::vcount(wfnb_igraph_net)-2)) / 2)
wfnb_betAll.norm <- (wfnb_betAll - min(wfnb_betAll))/(max(wfnb_betAll) - min(wfnb_betAll))
wfnb_dsAll <- igraph::similarity(method = "dice", wfnb_igraph_net, vids = igraph::V(wfnb_igraph_net), mode = "all")
wfnb_igraph_net <- igraph::set_vertex_attr(wfnb_igraph_net, "degree", index = igraph::V(wfnb_igraph_net), value = wfnb_degAll)
wfnb_igraph_net <- igraph::set_vertex_attr(wfnb_igraph_net, "betweenness", index = igraph::V(wfnb_igraph_net), value = wfnb_betAll.norm)
wfnb_igraph_net <- igraph::set_edge_attr(wfnb_igraph_net, "weight", index = igraph::E(wfnb_igraph_net), value = NULL)

nfwb_degAll <- igraph::degree(nfwb_igraph_net, v = igraph::V(nfwb_igraph_net), mode = "all")
nfwb_betAll <- igraph::betweenness(nfwb_igraph_net, v = igraph::V(nfwb_igraph_net), directed = FALSE) / (((igraph::vcount(nfwb_igraph_net) - 1) * (igraph::vcount(nfwb_igraph_net)-2)) / 2)
nfwb_betAll.norm <- (nfwb_betAll - min(nfwb_betAll))/(max(nfwb_betAll) - min(nfwb_betAll))
nfwb_dsAll <- igraph::similarity(method = "dice", nfwb_igraph_net, vids = igraph::V(nfwb_igraph_net), mode = "all")
nfwb_igraph_net <- igraph::set_vertex_attr(nfwb_igraph_net, "degree", index = igraph::V(nfwb_igraph_net), value = nfwb_degAll)
nfwb_igraph_net <- igraph::set_vertex_attr(nfwb_igraph_net, "betweenness", index = igraph::V(nfwb_igraph_net), value = nfwb_betAll.norm)
nfwb_igraph_net <- igraph::set_edge_attr(nfwb_igraph_net, "weight", index = igraph::E(nfwb_igraph_net), value = NULL)

# Create network in Cytoscape
createNetworkFromIgraph(wfb_igraph_net, new.title='Fusarium Brome')
createNetworkFromIgraph(nfb_igraph_net, new.title='No Fusarium No Brome')
createNetworkFromIgraph(wfnb_igraph_net, new.title='Fusarium No Brome')
createNetworkFromIgraph(nfwb_igraph_net, new.title='No Fusarium Brome')

# Import network tables from Cytoscape
wfb_at <- read.csv("Fusarium Brome default node.csv")
nfb_at <- read.csv("No Fusarium Brome default node.csv")
wfnb_at <- read.csv("Fusarium No Brome default node.csv")
nfwb_at <- read.csv("No Fusarium Brome default node.csv")

wfb_at$Treatment <- "wfb"
nfb_at$Treatment <- "nfb"
wfnb_at$Treatment <- "wfnb"
nfwb_at$Treatment <- "nfwb"

mdf <- rbind(wfb_at, nfb_at, wfnb_at, nfwb_at)

# Plots
library(tidyr)
library(RColorBrewer)
library(ggpubr)
library(phylosmith)

#df_melt <- mdf[,c(1:5,26,27)] # Subset to data to compare.
df_melt <- mdf[,c(3:5,27)] 
df_melt <- df_melt %>%
  mutate(Treatment = dplyr::recode(Treatment, nfb = "No F. pseudograminearum No B. tectorum", wfb = "F. pseudograminearum B. tectorum", 
                            nfwb = "No F. pseudograminearum B. tectorum", wfnb = "F. pseudograminearum No B. tectorum"))

df_melt$Treatment <- factor(df_melt$Treatment, levels = c("No F. pseudograminearum No B. tectorum", "No F. pseudograminearum B. tectorum",  "F. pseudograminearum B. tectorum", "F. pseudograminearum No B. tectorum"))

df_long <- df_melt %>% 
  pivot_longer(cols = 1:3, names_to = "Variable", values_to = "Value")

# Get mean value for each variable by treatment
dfx<-aggregate(df_long[,3],by=list(df_long$Treatment, df_long$Variable),FUN=mean)
(0.32-0.28)/0.32 # ClosenessCentrality
(9.5-7.4)/7.4 # Degree: 28% reduction

my_comp <- list(c("No F. pseudograminearum No B. tectorum", "F. pseudograminearum B. tectorum"), c("No F. pseudograminearum No B. tectorum", "No F. pseudograminearum B. tectorum"), c("No F. pseudograminearum No B. tectorum", "F. pseudograminearum No B. tectorum"))

v1 <- ggviolin(df_long, x = "Treatment", y = "Value", fill = "Treatment", orientation = "horizontal", palette = brewer.pal(4,"Dark2"),
               add = "boxplot", add.params = list(fill = "white")) + 
  stat_compare_means(comparisons = my_comp, label = "p.signif")+ # Add significance levels
  facet_wrap(~Variable, scales = "free") + theme(legend.position = "bottom", axis.text.y = element_blank()) +
  theme(plot.margin = margin(t = 0, r = 0.5, b = 0, l = 0.5, unit = "cm"), panel.spacing = unit(0.5, "lines")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
v1

pdf("Figures/Network Properties Fusarium Brome.pdf", width = 10, height = 6)
v1 
dev.off()

m16$sample_table$Fus_Brome <- paste0(m16$sample_table$F..pseudograminearumPresence, " ", m16$sample_table$B.tectorum)
m16$sample_table$Fus_Brome <- factor(m16$sample_table$Fus_Brome, levels = c( "F. pseudograminearum B. tectorum", "F. pseudograminearum No B. tectorum", "No F. pseudograminearum B. tectorum", "No F. pseudograminearum No B. tectorum"))
m16$tidy_dataset()
m16$sample_table$Fus_Brome

trt <- m16$merge_samples(group = "Fus_Brome")
trt$sample_table$SampleID <- factor(trt$sample_table$SampleID, levels = c( "F. pseudograminearum B. tectorum", "F. pseudograminearum No B. tectorum", "No F. pseudograminearum B. tectorum", "No F. pseudograminearum No B. tectoru"))
trt$tidy_dataset()

t3 <- trans_venn$new(trt, ratio = "numratio")

g1 <- t3$plot_bar(left_plot = TRUE, bottom_height = 0.5, left_width = 0.15, up_bar_fill = "grey50", left_bar_fill = "grey50", bottom_point_color = "black")
g1

pdf("Figures/Fusarium Brome Upset Plot.pdf", width = 10, height = 6)
g1
dev.off()


# Hub taxa calculated with Cytohubba in Cytoscape and .csvs imported into R for analysis
wfb_hub <- read.csv("Fusarium Brome_MCC_top20 default node.csv")
nfb_hub <- read.csv("No Fusarium No Brome_MCC_top20 default node.csv")
wfnb_hub <- read.csv("Fusarium No Brome_MCC_top20 default node.csv")
nfwb_hub <- read.csv("No Fusarium Brome_MCC_top20 default node.csv")

wfb_hub$Treatment <- "With Fusarium and Brome"
nfb_hub$Treatment <- "No Fusarium or Brome"
wfnb_hub$Treatment <- "With Fusarium No Brome"
nfwb_hub$Treatment <- "No Fusarium With Brome"

hub_df <- rbind(wfb_hub, nfb_hub, wfnb_hub, nfwb_hub)
write.csv(hub_df, "Hub Taxa Fusarium Brome.csv")

unique(hub_df$Phylum)
hub_df$Phylum <- factor(hub_df$Phylum, levels = c("Acidobacteria", "Actinobacteria",  "Bacteroidetes", "Firmicutes", "Armatimonadetes", "Proteobacteria", "Cyanobacteria"))
ggbarplot(hub_df, x = "Treatment", y= "RelativeAbundance", fill = "Phylum") + ylab("Relative Abundance (%)")

unique(hub_df$Family)
ggbarplot(hub_df, x = "Treatment", y= "RelativeAbundance", fill = "Family") + ylab("Relative Abundance (%)")

# Venn Diagram
library(ggvenn)

# Phylum
x2 <- list('With Fusarium and Brome' = wfb_hub$Phylum, 'No Fusarium or Brome' = nfb_hub$Phylum, "With Fusarium No Brome" = wfnb_hub$Phylum, 'No Fusarium With Brome' = nfwb_hub$Phylum)

pdf("Figures/Venn Phylum Hub.pdf", width = 10, height = 8)
ggvenn(
  x2, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#CD80C2FF" ,"#ED73C2FF"),
  stroke_size = 0, set_name_size = 4, show_outside = "always"
)
dev.off() 

# Family
x3 <- list('Fusarium Brome' = wfb_hub$Family, 'No Fusarium No Brome' = nfb_hub$Family, "Fusarium No Brome" = wfnb_hub$Family, 'No Fusarium Brome' = nfwb_hub$Family)

pdf("Figures/Venn Family Hub.pdf", width = 10, height = 8)
ggvenn(
  x3, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#CD80C2FF" ,"#ED73C2FF"),
  stroke_size = 0, set_name_size = 4, show_outside = "always"
)
dev.off() 

# Genus
x4 <- list('Fusarium Brome' = wfb_hub$Genus, 'No Fusarium No Brome' = nfb_hub$Genus, "Fusarium No Brome" = wfnb_hub$Genus, 'No Fusarium Brome' = nfwb_hub$Genus)

pdf("Figures/Venn Genus Hub.pdf", width = 10, height = 8)
ggvenn(
  x4, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#CD80C2FF" ,"#ED73C2FF"),
  stroke_size = 0, set_name_size = 4, show_outside = "always"
)
dev.off() 


# UpSet Plot of network nodes
ps <- readRDS("../../merged_phyloseq.RDS")
sd <- sample_data(ps)

udf <- rbind(wfb_at, nfb_at, wfnb_at, nfwb_at)
# Phylum level
udP <- udf[, c(18,27)]
# Family level
udF <- udf[, c(7,27)]

# Phylum level
psp <- tax_glom(ps, taxrank = "Phylum")
tax <- data.frame(tax_table(psp))
tax$ASV <- rownames(tax)

uudP <- unique(udP)
taxm <- merge(uudP[1], tax[1:8], by = "Phylum", all = TRUE)
taxm <- distinct(taxm)
taxm <- subset(taxm, ASV != "NA")
rownames(taxm) <- taxm$ASV
taxa <- taxm[,c(1,2)]
taxa <- taxa[,c(2,1)]
rows <- rownames(taxa)
taxa <- as.matrix(taxa)

otu <- as.matrix(otu_table(ps))
otu <- otu[rownames(otu) %in% rows]

sd$Fus_brom <- paste0(sd$FusariumPresence, " ", sd$Weeds)

psp_sub <- phyloseq(tax_table(taxa), otu_table(otu), sample_data(sd))
psp_sub

# Family level
psf <- tax_glom(ps, taxrank = "Family")
tax <- data.frame(tax_table(psf))
tax$ASV <- rownames(tax)

uudF <- unique(udF)
taxf <- merge(uudF[1], tax[1:8], by = "Family", all = TRUE)
taxf <- distinct(taxf)
taxf <- subset(taxf, ASV != "NA")
rownames(taxf) <- taxf$ASV
taxaf <- taxf[,c(1:2)]
taxaf <- taxaf[,c(2,1)]
rowsf <- rownames(taxaf)
taxaf <- as.matrix(taxaf)

otu <- as.matrix(otu_table(ps))
otu <- otu[rownames(otu) %in% rowsf]

psf_sub <- phyloseq(tax_table(taxaf), otu_table(otu), sample_data(sd))
psf_sub

mecoP <- phyloseq2meco(psp_sub)
mecoF <- phyloseq2meco(psf_sub)
mecoP$sample_table

mecoP$sample_table$Fus_brom <- paste0(mecoP$sample_table$F..pseudograminearumPresence, " ", mecoP$sample_table$B.tectorum)
mecoP$sample_table$Fus_brom <- factor(mecoP$sample_table$Fus_Brome, levels = c( "F. pseudograminearum B. tectorum", "F. pseudograminearum No B. tectorum", "No F. pseudograminearum B. tectorum", "No F. pseudograminearum No B. tectorum"))

mecoF$sample_table
mecoF$sample_table$Fus_brom <- paste0(mecoF$sample_table$F..pseudograminearumPresence, " ", mecoF$sample_table$B.tectorum)
mecoF$sample_table$Fus_brom <- factor(mecoF$sample_table$Fus_Brome, levels = c( "F. pseudograminearum B. tectorum", "F. pseudograminearum No B. tectorum", "No F. pseudograminearum B. tectorum", "No F. pseudograminearum No B. tectorum"))

mecoP$tidy_dataset()
mecoF$tidy_dataset()

tP <- mecoP$merge_samples("Fus_brom")
tF <- mecoF$merge_samples("Fus_brom")

tP$sample_table$SampleID <- factor(tP$sample_table$SampleID, levels = c( "F. pseudograminearum B. tectorum", "F. pseudograminearum No B. tectorum", "No F. pseudograminearum B. tectorum", "No F. pseudograminearum No B. tectorum"))
tP$tidy_dataset()

tvP <- trans_venn$new(dataset = tP)
tvF <- trans_venn$new(dataset = tF)

tvP$data_summary %<>% .[.[, 1] > 0, ]
tpP <- tvP$plot_bar(left_plot = TRUE, bottom_height = 0.5, left_width = 0.15, up_bar_fill = "grey50", left_bar_fill = "grey50", bottom_point_color = "black", bottom_y_text_size = 0) 
tpP 


tvF$data_summary %<>% .[.[, 1] > 0, ]
tpF <- tvF$plot_bar(left_plot = TRUE, bottom_height = 0.5, left_width = 0.15, up_bar_fill = "grey50", left_bar_fill = "grey50", bottom_point_color = "black")
tpF

# Convert UpSet plots to ggplots to arrange
library(ggplotify)
tgP <- as.ggplot(tpP)
tgP

np1 <- as.ggplot(np)

vp1 <- ggarrange(np1, tgP, nrow = 1, ncol = 2, labels = c("a", "b"), widths = c(1,0.5)) 
vp1

v11 <- v1 + theme(legend.title=element_blank())

vp2 <- ggarrange(v11, cp1, nrow = 1, ncol = 2, labels = c("c", "d"), widths = c(1,0.4), common.legend = TRUE, legend = "bottom")
vp2

ggarrange(vp1, vp2, nrow = 2, ncol = 1, heights = c(1.5, 1)) + theme(plot.margin = margin(0,0,0,0, "cm"))

pdf("Figures/Fig4.pdf", width = 12, height = 8)
ggarrange(vp1, vp2, nrow = 2, ncol = 1, heights = c(1.5, 1)) 
dev.off() 


# Combined Network Figure (Fig. S2)
library(png)
library(grid)
library(ggplot2)
library(ggpubr)

# import Gephi plots
nfb.img <- readPNG("Figures/nfb_net.png")
wfb.img <- readPNG("Figures/wfb_net.png")
wfnb.img <- readPNG("Figures/wfnb_net.png")
nfwb.img <- readPNG("Figures/nfwb_net.png")

nfbg <- rasterGrob(nfb.img, interpolate = TRUE)
wfbg <- rasterGrob(wfb.img, interpolate = TRUE)
wfnbg <- rasterGrob(wfnb.img, interpolate = TRUE)
nfwbg <- rasterGrob(nfwb.img, interpolate = TRUE)

p1 <- ggplot() + annotation_custom(nfbg) + theme_classic()
p2 <- ggplot() + annotation_custom(wfbg) + theme_classic()
p3 <- ggplot() + annotation_custom(wfnbg) + theme_classic()
p4 <- ggplot() + annotation_custom(nfwbg) + theme_classic()

ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"), heights = c(1,1))

pdf("Figures/Network Combined.pdf", width = 8, height = 8)
ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"), heights = c(1,1))
dev.off()


########## END Across Locations ##############

# By location year
bz22wf <- clone(m16)
bz22nf <- clone(m16)
bz23wf <- clone(m16)
bz23nf <- clone(m16)
sarc22wf <- clone(m16)
sarc22nf <- clone(m16)
sarc23wf <- clone(m16)
sarc23nf <- clone(m16)

bz22wf$sample_table %<>% subset(Location == "Bozeman" & Year == "2022" & FusariumPresence == "Fusarium")
bz22nf$sample_table %<>% subset(Location == "Bozeman" & Year == "2022" & FusariumPresence == "No Fusarium")
bz23wf$sample_table %<>% subset(Location == "Bozeman" & Year == "2023" & FusariumPresence == "Fusarium")
bz23nf$sample_table %<>% subset(Location == "Bozeman" & Year == "2023" & FusariumPresence == "No Fusarium")

sarc22wf$sample_table %<>% subset(Location == "Huntley" & Year == "2022" & FusariumPresence == "Fusarium")
sarc22nf$sample_table %<>% subset(Location == "Huntley" & Year == "2022" & FusariumPresence == "No Fusarium")
sarc23wf$sample_table %<>% subset(Location == "Huntley" & Year == "2023" & FusariumPresence == "Fusarium")
sarc23nf$sample_table %<>% subset(Location == "Huntley" & Year == "2023" & FusariumPresence == "No Fusarium")

bz22wf$tidy_dataset()
bz22nf$tidy_dataset()
bz23wf$tidy_dataset()
bz23nf$tidy_dataset()

sarc22wf$tidy_dataset()
sarc22nf$tidy_dataset()
sarc23wf$tidy_dataset()
sarc23nf$tidy_dataset()

# Filter to 50% prevalence
#bz22wf$filter_taxa(freq = 0.5)


bz22wf_net <- trans_network$new(dataset = bz22wf, cor_method = NULL, taxa_level = "OTU", filter_thres = 0.0002)
bz22wf_net$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb", icov.select.params(ncores = 40), 
                       add_taxa_name = c("Kingdom", "Phylum", "Family", "Genus"), verbose=TRUE)
saveRDS(bz22wf_net, "bz22wf_net.RDS")

bz22nf_net <- trans_network$new(dataset = bz22nf, cor_method = NULL, taxa_level = "OTU", filter_thres = 0.0002)
bz22nf_net$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb", icov.select.params(ncores = 40), 
                       add_taxa_name = c("Kingdom", "Phylum", "Family", "Genus"), verbose=TRUE)
saveRDS(bz22nf_net, "bz22nf_net.RDS")

bz23wf_net <- trans_network$new(dataset = bz23wf, cor_method = NULL, taxa_level = "OTU", filter_thres = 0.0004)
bz23wf_net$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb", icov.select.params(ncores = 40), 
                       add_taxa_name = c("Kingdom", "Phylum", "Family", "Genus"), verbose=TRUE)
saveRDS(bz23wf_net, "bz23wf_net.RDS")

bz23nf_net <- trans_network$new(dataset = bz23nf, cor_method = NULL, taxa_level = "OTU", filter_thres = 0.0004)
bz23nf_net$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb", icov.select.params(ncores = 40), 
                       add_taxa_name = c("Kingdom", "Phylum", "Family", "Genus"), verbose=TRUE)
saveRDS(bz23nf_net, "bz23nf_net.RDS")

sarc22wf_net <- trans_network$new(dataset = sarc22wf, cor_method = NULL, taxa_level = "OTU", filter_thres = 0.0002)
sarc22wf_net$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb", icov.select.params(ncores = 40), 
                         add_taxa_name = c("Kingdom", "Phylum", "Family", "Genus"), verbose=TRUE)
saveRDS(sarc22wf_net, "sarc22wf_net.RDS")

sarc22nf_net <- trans_network$new(dataset = sarc22nf, cor_method = NULL, taxa_level = "OTU", filter_thres = 0.0003)
sarc22nf_net$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb", icov.select.params(ncores = 40), 
                         add_taxa_name = c("Kingdom", "Phylum", "Family", "Genus"), verbose=TRUE)
saveRDS(sarc22nf_net, "sarc22nf_net.RDS")

sarc23wf_net <- trans_network$new(dataset = sarc23wf, cor_method = NULL, taxa_level = "OTU", filter_thres = 0.00005)
sarc23wf_net$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb", icov.select.params(ncores = 40), 
                         add_taxa_name = c("Kingdom", "Phylum", "Family", "Genus"), verbose=TRUE)
saveRDS(sarc23wf_net, "sarc23wf_net.RDS")

sarc23nf_net <- trans_network$new(dataset = sarc23nf, cor_method = NULL, taxa_level = "OTU", filter_thres = 0.00005)
sarc23nf_net$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb", icov.select.params(ncores = 40), 
                         add_taxa_name = c("Kingdom", "Phylum", "Family", "Genus"), verbose=TRUE)
saveRDS(sarc23nf_net, "sarc23nf_net.RDS")

# Import previously constructed networks
bz22wf_net <- readRDS("bz22wf_net.RDS")
bz22nf_net <- readRDS("bz22nf_net.RDS")
bz23wf_net <- readRDS("bz23wf_net.RDS")
bz23nf_net <- readRDS("bz23nf_net.RDS")
sarc22wf_net <- readRDS("sarc22wf_net.RDS")
sarc22nf_net <- readRDS("sarc22nf_net.RDS")
sarc23wf_net <- readRDS("sarc23wf_net.RDS")
sarc23nf_net <- readRDS("sarc23nf_net.RDS")

bz22wf_net$cal_module(method = "cluster_fast_greedy")
bz22nf_net$cal_module(method = "cluster_fast_greedy")
bz23wf_net$cal_module(method = "cluster_fast_greedy")
bz23nf_net$cal_module(method = "cluster_fast_greedy")
sarc22wf_net$cal_module(method = "cluster_fast_greedy")
sarc22nf_net$cal_module(method = "cluster_fast_greedy")
sarc23wf_net$cal_module(method = "cluster_fast_greedy")
sarc23nf_net$cal_module(method = "cluster_fast_greedy")

bz22wf_net$cal_network_attr()
bz22nf_net$cal_network_attr()
bz23wf_net$cal_network_attr()
bz23nf_net$cal_network_attr()
sarc22wf_net$cal_network_attr()
sarc22nf_net$cal_network_attr()
sarc23wf_net$cal_network_attr()
sarc23nf_net$cal_network_attr()

bz22wf_net$get_edge_table()
bz22nf_net$get_edge_table()
bz23wf_net$get_edge_table()
bz23nf_net$get_edge_table()
sarc22wf_net$get_edge_table()
sarc22nf_net$get_edge_table()
sarc23wf_net$get_edge_table()
sarc23nf_net$get_edge_table()

bz22wf_net$get_adjacency_matrix()
bz22nf_net$get_adjacency_matrix()
bz23wf_net$get_adjacency_matrix()
bz23nf_net$get_adjacency_matrix()
sarc22wf_net$get_adjacency_matrix()
sarc22nf_net$get_adjacency_matrix()
sarc23wf_net$get_adjacency_matrix()
sarc23nf_net$get_adjacency_matrix()

bz22wf_net$get_node_table(node_roles = TRUE)
bz22nf_net$get_node_table(node_roles = TRUE)
bz23wf_net$get_node_table(node_roles = TRUE)
bz23nf_net$get_node_table(node_roles = TRUE)
sarc22wf_net$get_node_table(node_roles = TRUE)
sarc22nf_net$get_node_table(node_roles = TRUE)
sarc23wf_net$get_node_table(node_roles = TRUE)
sarc23nf_net$get_node_table(node_roles = TRUE)

bz22wf_df <- bz22wf_net$res_network_attr
bz22nf_df <- bz22nf_net$res_network_attr
bz23wf_df <- bz23wf_net$res_network_attr
bz23nf_df <- bz23nf_net$res_network_attr
sarc22wf_df <- sarc22wf_net$res_network_attr
sarc22nf_df <- sarc22nf_net$res_network_attr
sarc23wf_df <- sarc23wf_net$res_network_attr
sarc23nf_df <- sarc23nf_net$res_network_attr

net_df <- cbind(bz22wf_df,bz22nf_df,bz23wf_df,bz23nf_df,sarc22wf_df,sarc22nf_df,sarc23wf_df,sarc23nf_df)
write.csv(net_df, "Network Attributes.csv")

bz22wf_net$plot_taxa_roles(use_type = 1)
bz22nf_net$plot_taxa_roles(use_type = 1)
bz23wf_net$plot_taxa_roles(use_type = 1)
bz23nf_net$plot_taxa_roles(use_type = 1)
sarc22wf_net$plot_taxa_roles(use_type = 1)
sarc22nf_net$plot_taxa_roles(use_type = 1)
sarc23wf_net$plot_taxa_roles(use_type = 1)
sarc23nf_net$plot_taxa_roles(use_type = 1)

bz22wf_net$plot_taxa_roles(use_type = 2)
bz22nf_net$plot_taxa_roles(use_type = 2)
bz23wf_net$plot_taxa_roles(use_type = 2)
bz23nf_net$plot_taxa_roles(use_type = 2)
sarc22wf_net$plot_taxa_roles(use_type = 2)
sarc22nf_net$plot_taxa_roles(use_type = 2)
sarc23wf_net$plot_taxa_roles(use_type = 2)
sarc23nf_net$plot_taxa_roles(use_type = 2)

bz22wf_net$cal_sum_links()
bz22nf_net$cal_sum_links()
bz23wf_net$cal_sum_links()
bz23nf_net$cal_sum_links()
sarc22wf_net$cal_sum_links()
sarc22nf_net$cal_sum_links()
sarc23wf_net$cal_sum_links()
sarc23nf_net$cal_sum_links()

bz22wf_net$plot_sum_links(method = "circlize", transparency = 0.2, plot_num = 6, annotationTrackHeight = circlize::mm_h(c(5, 5)))
bz22nf_net$plot_sum_links(method = "circlize", transparency = 0.2, plot_num = 6, annotationTrackHeight = circlize::mm_h(c(5, 5)))
bz23wf_net$plot_sum_links(method = "circlize", transparency = 0.2, plot_num = 6, annotationTrackHeight = circlize::mm_h(c(5, 5)))
bz23nf_net$plot_sum_links(method = "circlize", transparency = 0.2, plot_num = 6, annotationTrackHeight = circlize::mm_h(c(5, 5)))
sarc22wf_net$plot_sum_links(method = "circlize", transparency = 0.2, plot_num = 6, annotationTrackHeight = circlize::mm_h(c(5, 5)))
sarc22nf_net$plot_sum_links(method = "circlize", transparency = 0.2, plot_num = 6, annotationTrackHeight = circlize::mm_h(c(5, 5)))
sarc23wf_net$plot_sum_links(method = "circlize", transparency = 0.2, plot_num = 6, annotationTrackHeight = circlize::mm_h(c(5, 5)))
sarc23nf_net$plot_sum_links(method = "circlize", transparency = 0.2, plot_num = 6, annotationTrackHeight = circlize::mm_h(c(5, 5)))


bz22wf_net$plot_network(method = "ggraph", node_color="module")
bz22nf_net$plot_network(method = "ggraph", node_color="module")
bz23wf_net$plot_network(method = "ggraph", node_color="module")
bz23nf_net$plot_network(method = "ggraph", node_color="module")
sarc22wf_net$plot_network(method = "ggraph", node_color="module")
sarc22nf_net$plot_network(method = "ggraph", node_color="module")
sarc23wf_net$plot_network(method = "ggraph", node_color="module")
sarc23nf_net$plot_network(method = "ggraph", node_color="module")

bz22wf_net$save_network(filepath = "bz22wf_net.gexf")
bz22nf_net$save_network(filepath = "bz22nf_net.gexf")
bz23wf_net$save_network(filepath = "bz23wf_net.gexf")
bz23nf_net$save_network(filepath = "bz23nf_net.gexf")
sarc22wf_net$save_network(filepath = "sarc22wf_net.gexf")
sarc22nf_net$save_network(filepath = "sarc22nf_net.gexf")
sarc23wf_net$save_network(filepath = "sarc23wf_net.gexf")
sarc23nf_net$save_network(filepath = "sarc23nf_net.gexf")

# View directly in Cytoscape
# https://cytoscape.org/cytoscape-automation/for-scripters/R/notebooks/Network-functions-and-visualization.nb.html
library(RCy3)

bz22wf_igraph_net <- bz22wf_net$res_network
bz22nf_igraph_net <- bz22nf_net$res_network
bz23wf_igraph_net <- bz23wf_net$res_network
bz23nf_igraph_net <- bz23nf_net$res_network
sarc22wf_igraph_net <- sarc22wf_net$res_network
sarc22nf_igraph_net <- sarc22nf_net$res_network
sarc23wf_igraph_net <- sarc23wf_net$res_network
sarc23nf_igraph_net <- sarc23nf_net$res_network

bz22wf_degAll <- igraph::degree(bz22wf_igraph_net, v = igraph::V(bz22wf_igraph_net), mode = "all")
bz22nf_degAll <- igraph::degree(bz22nf_igraph_net, v = igraph::V(bz22nf_igraph_net), mode = "all")
bz23wf_degAll <- igraph::degree(bz23wf_igraph_net, v = igraph::V(bz23wf_igraph_net), mode = "all")
bz23nf_degAll <- igraph::degree(bz23nf_igraph_net, v = igraph::V(bz23nf_igraph_net), mode = "all")
sarc22wf_degAll <- igraph::degree(sarc22wf_igraph_net, v = igraph::V(sarc22wf_igraph_net), mode = "all")
sarc22nf_degAll <- igraph::degree(sarc22nf_igraph_net, v = igraph::V(sarc22nf_igraph_net), mode = "all")
sarc23wf_degAll <- igraph::degree(sarc23wf_igraph_net, v = igraph::V(sarc23wf_igraph_net), mode = "all")
sarc23nf_degAll <- igraph::degree(sarc23nf_igraph_net, v = igraph::V(sarc23nf_igraph_net), mode = "all")

bz22wf_betAll <- igraph::betweenness(bz22wf_igraph_net, v = igraph::V(bz22wf_igraph_net), directed = FALSE) / (((igraph::vcount(bz22wf_igraph_net) - 1) * (igraph::vcount(bz22wf_igraph_net)-2)) / 2)
bz22nf_betAll <- igraph::betweenness(bz22nf_igraph_net, v = igraph::V(bz22nf_igraph_net), directed = FALSE) / (((igraph::vcount(bz22nf_igraph_net) - 1) * (igraph::vcount(bz22nf_igraph_net)-2)) / 2)
bz23wf_betAll <- igraph::betweenness(bz23wf_igraph_net, v = igraph::V(bz23wf_igraph_net), directed = FALSE) / (((igraph::vcount(bz23wf_igraph_net) - 1) * (igraph::vcount(bz23wf_igraph_net)-2)) / 2)
bz23nf_betAll <- igraph::betweenness(bz23nf_igraph_net, v = igraph::V(bz23nf_igraph_net), directed = FALSE) / (((igraph::vcount(bz23nf_igraph_net) - 1) * (igraph::vcount(bz23nf_igraph_net)-2)) / 2)
sarc22wf_betAll <- igraph::betweenness(sarc22wf_igraph_net, v = igraph::V(sarc22wf_igraph_net), directed = FALSE) / (((igraph::vcount(sarc22wf_igraph_net) - 1) * (igraph::vcount(sarc22wf_igraph_net)-2)) / 2)
sarc22nf_betAll <- igraph::betweenness(sarc22nf_igraph_net, v = igraph::V(sarc22nf_igraph_net), directed = FALSE) / (((igraph::vcount(sarc22nf_igraph_net) - 1) * (igraph::vcount(sarc22nf_igraph_net)-2)) / 2)
sarc23wf_betAll <- igraph::betweenness(sarc23wf_igraph_net, v = igraph::V(sarc23wf_igraph_net), directed = FALSE) / (((igraph::vcount(sarc23wf_igraph_net) - 1) * (igraph::vcount(sarc23wf_igraph_net)-2)) / 2)
sarc23nf_betAll <- igraph::betweenness(sarc23nf_igraph_net, v = igraph::V(sarc23nf_igraph_net), directed = FALSE) / (((igraph::vcount(sarc23nf_igraph_net) - 1) * (igraph::vcount(sarc23nf_igraph_net)-2)) / 2)

bz22wf_betAll.norm <- (bz22wf_betAll - min(bz22wf_betAll))/(max(bz22wf_betAll) - min(bz22wf_betAll))
bz22nf_betAll.norm <- (bz22nf_betAll - min(bz22nf_betAll))/(max(bz22nf_betAll) - min(bz22nf_betAll))
bz23wf_betAll.norm <- (bz23wf_betAll - min(bz23wf_betAll))/(max(bz23wf_betAll) - min(bz23wf_betAll))
bz23nf_betAll.norm <- (bz23nf_betAll - min(bz23nf_betAll))/(max(bz23nf_betAll) - min(bz23nf_betAll))
sarc22wf_betAll.norm <- (sarc22wf_betAll - min(sarc22wf_betAll))/(max(sarc22wf_betAll) - min(sarc22wf_betAll))
sarc22nf_betAll.norm <- (sarc22nf_betAll - min(sarc22nf_betAll))/(max(sarc22nf_betAll) - min(sarc22nf_betAll))
sarc23wf_betAll.norm <- (sarc23wf_betAll - min(sarc23wf_betAll))/(max(sarc23wf_betAll) - min(sarc23wf_betAll))
sarc23nf_betAll.norm <- (sarc23nf_betAll - min(sarc23nf_betAll))/(max(sarc23nf_betAll) - min(sarc23nf_betAll))

bz22wf_dsAll <- igraph::similarity(method = "dice", bz22wf_igraph_net, vids = igraph::V(bz22wf_igraph_net), mode = "all")
bz22nf_dsAll <- igraph::similarity(method = "dice", bz22nf_igraph_net, vids = igraph::V(bz22nf_igraph_net), mode = "all")
bz23wf_dsAll <- igraph::similarity(method = "dice", bz23wf_igraph_net, vids = igraph::V(bz23wf_igraph_net), mode = "all")
bz23nf_dsAll <- igraph::similarity(method = "dice", bz23nf_igraph_net, vids = igraph::V(bz23nf_igraph_net), mode = "all")
sarc22wf_dsAll <- igraph::similarity(method = "dice", sarc22wf_igraph_net, vids = igraph::V(sarc22wf_igraph_net), mode = "all")
sarc22nf_dsAll <- igraph::similarity(method = "dice", sarc22nf_igraph_net, vids = igraph::V(sarc22nf_igraph_net), mode = "all")
sarc23wf_dsAll <- igraph::similarity(method = "dice", sarc23wf_igraph_net, vids = igraph::V(sarc23wf_igraph_net), mode = "all")
sarc23nf_dsAll <- igraph::similarity(method = "dice", sarc23nf_igraph_net, vids = igraph::V(sarc23nf_igraph_net), mode = "all")

bz22wf_igraph_net <- igraph::set_vertex_attr(bz22wf_igraph_net, "degree", index = igraph::V(bz22wf_igraph_net), value = bz22wf_degAll)
bz22nf_igraph_net <- igraph::set_vertex_attr(bz22nf_igraph_net, "degree", index = igraph::V(bz22nf_igraph_net), value = bz22nf_degAll)
bz23wf_igraph_net <- igraph::set_vertex_attr(bz23wf_igraph_net, "degree", index = igraph::V(bz23wf_igraph_net), value = bz23wf_degAll)
bz23nf_igraph_net <- igraph::set_vertex_attr(bz23nf_igraph_net, "degree", index = igraph::V(bz23nf_igraph_net), value = bz23nf_degAll)
sarc22wf_igraph_net <- igraph::set_vertex_attr(sarc22wf_igraph_net, "degree", index = igraph::V(sarc22wf_igraph_net), value = sarc22wf_degAll)
sarc22nf_igraph_net <- igraph::set_vertex_attr(sarc22nf_igraph_net, "degree", index = igraph::V(sarc22nf_igraph_net), value = sarc22nf_degAll)
sarc23wf_igraph_net <- igraph::set_vertex_attr(sarc23wf_igraph_net, "degree", index = igraph::V(sarc23wf_igraph_net), value = sarc23wf_degAll)
sarc23nf_igraph_net <- igraph::set_vertex_attr(sarc23nf_igraph_net, "degree", index = igraph::V(sarc23nf_igraph_net), value = sarc23nf_degAll)

bz22wf_igraph_net <- igraph::set_vertex_attr(bz22wf_igraph_net, "betweenness", index = igraph::V(bz22wf_igraph_net), value = bz22wf_betAll.norm)
bz22nf_igraph_net <- igraph::set_vertex_attr(bz22nf_igraph_net, "betweenness", index = igraph::V(bz22nf_igraph_net), value = bz22nf_betAll.norm)
bz23wf_igraph_net <- igraph::set_vertex_attr(bz23wf_igraph_net, "betweenness", index = igraph::V(bz23wf_igraph_net), value = bz23wf_betAll.norm)
bz23nf_igraph_net <- igraph::set_vertex_attr(bz23nf_igraph_net, "betweenness", index = igraph::V(bz23nf_igraph_net), value = bz23nf_betAll.norm)
sarc22wf_igraph_net <- igraph::set_vertex_attr(sarc22wf_igraph_net, "betweenness", index = igraph::V(sarc22wf_igraph_net), value = sarc22wf_betAll.norm)
sarc22nf_igraph_net <- igraph::set_vertex_attr(sarc22nf_igraph_net, "betweenness", index = igraph::V(sarc22nf_igraph_net), value = sarc22nf_betAll.norm)
sarc23wf_igraph_net <- igraph::set_vertex_attr(sarc23wf_igraph_net, "betweenness", index = igraph::V(sarc23wf_igraph_net), value = sarc23wf_betAll.norm)
sarc23nf_igraph_net <- igraph::set_vertex_attr(sarc23nf_igraph_net, "betweenness", index = igraph::V(sarc23nf_igraph_net), value = sarc23nf_betAll.norm)

bz22wf_igraph_net <- igraph::set_edge_attr(bz22wf_igraph_net, "weight", index = igraph::E(bz22wf_igraph_net), value = NULL)
bz22nf_igraph_net <- igraph::set_edge_attr(bz22nf_igraph_net, "weight", index = igraph::E(bz22nf_igraph_net), value = NULL)
bz23wf_igraph_net <- igraph::set_edge_attr(bz23wf_igraph_net, "weight", index = igraph::E(bz23wf_igraph_net), value = NULL)
bz23nf_igraph_net <- igraph::set_edge_attr(bz23nf_igraph_net, "weight", index = igraph::E(bz23nf_igraph_net), value = NULL)
sarc22wf_igraph_net <- igraph::set_edge_attr(sarc22wf_igraph_net, "weight", index = igraph::E(sarc22wf_igraph_net), value = NULL)
sarc22nf_igraph_net <- igraph::set_edge_attr(sarc22nf_igraph_net, "weight", index = igraph::E(sarc22nf_igraph_net), value = NULL)
sarc23wf_igraph_net <- igraph::set_edge_attr(sarc23wf_igraph_net, "weight", index = igraph::E(sarc23wf_igraph_net), value = NULL)
sarc23nf_igraph_net <- igraph::set_edge_attr(sarc23nf_igraph_net, "weight", index = igraph::E(sarc23nf_igraph_net), value = NULL)

# BZ network in Cytoscape
createNetworkFromIgraph(bz22wf_igraph_net,new.title='Bozeman 2022 Fusarium')
createNetworkFromIgraph(bz22nf_igraph_net,new.title='Bozeman 2022 No Fusarium')
createNetworkFromIgraph(bz23wf_igraph_net,new.title='Bozeman 2023 Fusarium')
createNetworkFromIgraph(bz23nf_igraph_net,new.title='Bozeman 2023 No Fusarium')
createNetworkFromIgraph(sarc22wf_igraph_net,new.title='Huntley 2022 Fusarium')
createNetworkFromIgraph(sarc22nf_igraph_net,new.title='Huntley 2022 No Fusarium')
createNetworkFromIgraph(sarc23wf_igraph_net,new.title='Huntley 2023 Fusarium')
createNetworkFromIgraph(sarc23nf_igraph_net,new.title='Huntley 2023 No Fusarium')

# The following style settings will be applied to all graphs
setNodeShapeDefault(new.shape = "ELLIPSE", style.name = "default")
lockNodeDimensions(TRUE)
setNodeSizeMapping('RelativeAbundance', c(min(bz22wf_betAll.norm), mean(bz22wf_betAll.norm), max(bz22wf_betAll.norm)), c(1, 10, 100))

# Import network tables from Cytoscape
bz22wf_at <- read.csv("Bozeman 2022 Fusarium default node.csv")
bz22nf_at <- read.csv("Bozeman 2022 No Fusarium default node.csv")
bz23wf_at <- read.csv("Bozeman 2023 Fusarium default node.csv")
bz23nf_at <- read.csv("Bozeman 2023 No Fusarium default node.csv")
sarc22wf_at <- read.csv("Huntley 2022 Fusarium default node.csv")
sarc22nf_at <- read.csv("Huntley 2022 No Fusarium default node.csv")
sarc23wf_at <- read.csv("Huntley 2023 Fusarium default node.csv")
sarc23nf_at <- read.csv("Huntley 2023 No Fusarium default node.csv")

bz22wf_at$Treatment <- "bz22wf"
bz22nf_at$Treatment <- "bz22nf"
bz23wf_at$Treatment <- "bz23wf"
bz23nf_at$Treatment <- "bz23nf"
sarc22wf_at$Treatment <- "sarc22wf"
sarc22nf_at$Treatment <- "sarc22nf"
sarc23wf_at$Treatment <- "sarc23wf"
sarc23nf_at$Treatment <- "sarc23nf"

mdf <- rbind(bz22wf_at, bz22nf_at, bz23wf_at, bz23nf_at, sarc22wf_at, sarc22nf_at, sarc23wf_at, sarc23nf_at)

# Plots
library(tidyr)
library(RColorBrewer)
library(ggpubr)
library(phylosmith)

df_melt <- mdf[,c(1,3:5,26,27)] # Subset to data to compare.

df_long <- df_melt %>% 
  pivot_longer(cols = 1:5, names_to = "Variable", values_to = "Value")

my_comp <- list(c("bz22wf_at", "bz22nf_at", "bz23wf_at", "bz23nf_at", "sarc22wf_at", "sarc22nf_at", "sarc23wf_at", "sarc23nf_at"))
v1 <- ggviolin(df_long, x = "Treatment", y = "Value", fill = "Treatment", orientation = "horizontal", 
               add = "boxplot", add.params = list(fill = "white")) + 
  stat_compare_means(comparisons = my_comp, label = "p.signif")+ # Add significance levels
  facet_wrap(~Variable, scales = "free") + theme(legend.position = "none") +
  theme(plot.margin = margin(t = 0, r = 0.5, b = 0, l = 0, unit = "cm")) 
v1 



# By location
m16 <- readRDS("../m16S_combined.RDS")
sd <- as.data.frame(m16$sample_table)

bz_wf <- clone(m16)
bz_nf <- clone(m16)
sarc_wf <- clone(m16)
sarc_nf <- clone(m16)

bz_wf$sample_table %<>% subset(Location == "Bozeman" & FusariumPresence == "Fusarium")
bz_nf$sample_table %<>% subset(Location == "Bozeman" & FusariumPresence == "No Fusarium")
sarc_wf$sample_table %<>% subset(Location == "Huntley" & FusariumPresence == "Fusarium")
sarc_nf$sample_table %<>% subset(Location == "Huntley" & FusariumPresence == "No Fusarium")

bz_wf$tidy_dataset()
bz_nf$tidy_dataset()
sarc_wf$tidy_dataset()
sarc_nf$tidy_dataset()

bzwf_net <- trans_network$new(dataset = bz_wf, cor_method = NULL, taxa_level = "OTU", filter_thres = 0.0005)
bzwf_net$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb", icov.select.params(ncores = 40), 
                       add_taxa_name = c("Kingdom", "Phylum", "Family", "Genus"), verbose=TRUE)
saveRDS(bzwf_net, "bzwf_net.RDS")

bznf_net <- trans_network$new(dataset = bz_nf, cor_method = NULL, taxa_level = "OTU", filter_thres = 0.0005)
bznf_net$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb", icov.select.params(ncores = 40), 
                     add_taxa_name = c("Kingdom", "Phylum", "Family", "Genus"), verbose=TRUE)
saveRDS(bznf_net, "bznf_net.RDS")

sarcwf_net <- trans_network$new(dataset = sarc_wf, cor_method = NULL, taxa_level = "OTU", filter_thres = 0.0005)
sarcwf_net$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb", icov.select.params(ncores = 40), 
                     add_taxa_name = c("Kingdom", "Phylum", "Family", "Genus"), verbose=TRUE)
saveRDS(sarcwf_net, "sarcwf_net.RDS")

sarcnf_net <- trans_network$new(dataset = sarc_nf, cor_method = NULL, taxa_level = "OTU", filter_thres = 0.0005)
sarcnf_net$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb", icov.select.params(ncores = 40), 
                       add_taxa_name = c("Kingdom", "Phylum", "Family", "Genus"), verbose=TRUE)
saveRDS(sarcnf_net, "sarcnf_net.RDS")

bzwf_net <- readRDS("bzwf_net.RDS")
bznf_net <- readRDS("bznf_net.RDS")
sarcwf_net <- readRDS("sarcwf_net.RDS")
sarcnf_net <- readRDS("sarcnf_net.RDS")

bzwf_net$cal_module(method = "cluster_fast_greedy")
bznf_net$cal_module(method = "cluster_fast_greedy")
sarcwf_net$cal_module(method = "cluster_fast_greedy")
sarcnf_net$cal_module(method = "cluster_fast_greedy")

bzwf_net$cal_network_attr()
bznf_net$cal_network_attr()
sarcwf_net$cal_network_attr()
sarcnf_net$cal_network_attr()

bzwf_net$save_network(filepath = "bzwf_net.gexf")
bznf_net$save_network(filepath = "bznf_net.gexf")
sarcwf_net$save_network(filepath = "sarcwf_net.gexf")
sarcnf_net$save_network(filepath = "sarcnf_net.gexf")

bzwf_net$res_network_attr
bznf_net$res_network_attr
sarcwf_net$res_network_attr
sarcnf_net$res_network_attr

bzwf_net$cal_eigen()
bznf_net$cal_eigen()
sarcwf_net$cal_eigen()
sarcnf_net$cal_eigen()

dim(bz_wf$tax_table)
bz_wf1 <- bz_wf$filter_taxa(rel_abund = 0.001, freq = 0.1)
dim(bz_wf1$tax_table)

bzwf_env <- trans_env$new(dataset = bz_wf, add_data = bz_wf$sample_table[, c(12:16,19)])
bznf_env <- trans_env$new(dataset = bz_nf, add_data = bz_nf$sample_table[, c(12:16,19)])
sarcwf_env <- trans_env$new(dataset = sarc_wf, add_data = sarc_wf$sample_table[, c(12:16,19)])
sarcnf_env <- trans_env$new(dataset = sarc_nf, add_data = sarc_nf$sample_table[, c(12:16,19)])

bzwf_env$cal_cor(add_abund_table = bzwf_net$res_eigen)
bznf_env$cal_cor(add_abund_table = bznf_net$res_eigen)
sarcwf_env$cal_cor(add_abund_table = sarcwf_net$res_eigen)
sarcnf_env$cal_cor(add_abund_table = sarcnf_net$res_eigen)

bzwf_env$plot_cor(color_palette = (RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")), text_y_order = c("M1", "M2", "M3","M4","M5","M6","M7","M8","M9","M10","M11","M12","M13","M14","M15","M16","M17","M18"))
bznf_env$plot_cor(color_palette = (RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")), text_y_order = c("M1", "M2", "M3","M4","M5","M6","M7","M8","M9","M10","M11","M12","M13","M14"))
sarcwf_env$plot_cor(color_palette = (RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")), text_y_order = c("M1", "M2", "M3","M4","M5","M6","M7","M8","M9","M10","M11","M12","M13","M14","M15"))
sarcnf_env$plot_cor(color_palette = (RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")), text_y_order = c("M1", "M2", "M3","M4","M5","M6","M7","M8","M9","M10","M11","M12","M13","M14","M15","M16","M17"))



bzwf_env$plot_cor()



# View directly in Cytoscape
# https://cytoscape.org/cytoscape-automation/for-scripters/R/notebooks/Network-functions-and-visualization.nb.html
library(RCy3)

bzwf_igraph_net <- bzwf_net$res_network
bznf_igraph_net <- bznf_net$res_network
sarcwf_igraph_net <- sarcwf_net$res_network
sarcnf_igraph_net <- sarcnf_net$res_network

bzwf_degAll <- igraph::degree(bzwf_igraph_net, v = igraph::V(bzwf_igraph_net), mode = "all")
bznf_degAll <- igraph::degree(bznf_igraph_net, v = igraph::V(bznf_igraph_net), mode = "all")
sarcwf_degAll <- igraph::degree(sarcwf_igraph_net, v = igraph::V(sarcwf_igraph_net), mode = "all")
sarcnf_degAll <- igraph::degree(sarcnf_igraph_net, v = igraph::V(sarcnf_igraph_net), mode = "all")

bzwf_betAll <- igraph::betweenness(bzwf_igraph_net, v = igraph::V(bzwf_igraph_net), directed = FALSE) / (((igraph::vcount(bzwf_igraph_net) - 1) * (igraph::vcount(bzwf_igraph_net)-2)) / 2)
bznf_betAll <- igraph::betweenness(bznf_igraph_net, v = igraph::V(bznf_igraph_net), directed = FALSE) / (((igraph::vcount(bznf_igraph_net) - 1) * (igraph::vcount(bznf_igraph_net)-2)) / 2)
sarcwf_betAll <- igraph::betweenness(sarcwf_igraph_net, v = igraph::V(sarcwf_igraph_net), directed = FALSE) / (((igraph::vcount(sarcwf_igraph_net) - 1) * (igraph::vcount(sarcwf_igraph_net)-2)) / 2)
sarcnf_betAll <- igraph::betweenness(sarcnf_igraph_net, v = igraph::V(sarcnf_igraph_net), directed = FALSE) / (((igraph::vcount(sarcnf_igraph_net) - 1) * (igraph::vcount(sarcnf_igraph_net)-2)) / 2)

bzwf_betAll.norm <- (bzwf_betAll - min(bzwf_betAll))/(max(bzwf_betAll) - min(bzwf_betAll))
bznf_betAll.norm <- (bznf_betAll - min(bznf_betAll))/(max(bznf_betAll) - min(bznf_betAll))
sarcwf_betAll.norm <- (sarcwf_betAll - min(sarcwf_betAll))/(max(sarcwf_betAll) - min(sarcwf_betAll))
sarcnf_betAll.norm <- (sarcnf_betAll - min(sarcnf_betAll))/(max(sarcnf_betAll) - min(sarcnf_betAll))

bzwf_dsAll <- igraph::similarity(method = "dice", bzwf_igraph_net, vids = igraph::V(bzwf_igraph_net), mode = "all")
bznf_dsAll <- igraph::similarity(method = "dice", bznf_igraph_net, vids = igraph::V(bznf_igraph_net), mode = "all")
sarcwf_dsAll <- igraph::similarity(method = "dice", sarcwf_igraph_net, vids = igraph::V(sarcwf_igraph_net), mode = "all")
sarcnf_dsAll <- igraph::similarity(method = "dice", sarcnf_igraph_net, vids = igraph::V(sarcnf_igraph_net), mode = "all")


# BZ network in Cytoscape
createNetworkFromIgraph(bzwf_igraph_net,new.title='Bozeman Fusarium')
createNetworkFromIgraph(bznf_igraph_net,new.title='Bozeman No Fusarium')
createNetworkFromIgraph(sarcwf_igraph_net,new.title='Huntley Fusarium')
createNetworkFromIgraph(sarcnf_igraph_net,new.title='Huntley No Fusarium')

# Import network tables from Cytoscape
bzwf_at <- read.csv("Bozeman Fusarium default node.csv")
bznf_at <- read.csv("Bozeman No Fusarium default node.csv")
sarcwf_at <- read.csv("Huntley Fusarium default node.csv")
sarcnf_at <- read.csv("Huntley No Fusarium default node.csv")

bzwf_at$Treatment <- "Bozeman Fusarium"
bznf_at$Treatment <- "Bozeman No Fusarium"
sarcwf_at$Treatment <- "Huntley Fusarium"
sarcnf_at$Treatment <- "Huntley No Fusarium"

mdf <- rbind(bzwf_at, bznf_at, sarcwf_at, sarcnf_at)

#df_melt <- mdf[,c(1,3:5,26,27)] # Subset to data to compare.
df_melt <- mdf[,c(1:5,26,27)]
df_melt <- subset(df_melt, ClosenessCentrality < 0.75)
df_melt <- subset(df_melt, ClusteringCoefficient < 0.8)

df_long <- df_melt %>% 
  pivot_longer(cols = 1:6, names_to = "Variable", values_to = "Value")

my_comp <- list(c("Bozeman Fusarium", "Bozeman No Fusarium"), c("Huntley Fusarium", "Huntley No Fusarium"))
v1 <- ggviolin(df_long, x = "Treatment", y = "Value", fill = "Treatment", ylab = FALSE, xlab = FALSE, orientation = "horizontal", 
               add = "boxplot", add.params = list(fill = "white")) + 
  facet_wrap(~Variable, scales = "free") + theme(legend.position = "bottom") +
  stat_compare_means(comparisons = my_comp, label = "p.signif")+ # Add significance levels
  theme(plot.margin = margin(t = 0, r = 0.5, b = 0, l = 0, unit = "cm")) 
v1 + scale_x_discrete(labels = NULL) 

pdf("Figures/Network Properties Fusarium.pdf", width = 10, height = 8)
v1 + scale_x_discrete(labels = NULL)
dev.off()

# Module comparison
m1 <- subset(mdf, module == "M1")
m2 <- subset(mdf, module == "M2")
m3 <- subset(mdf, module == "M3")
m4 <- subset(mdf, module == "M4")

comp <- list(c("Acidobacteria", "Actinobacteria"), c("Acidobacteria","Planctomycetes"))
ggboxplot(m1, x="Treatment", y="RelativeAbundance", fill = "Phylum")
ggboxplot(m2, x="Treatment", y="RelativeAbundance", fill = "Phylum")
ggboxplot(m3, x="Treatment", y="RelativeAbundance", fill = "Phylum") + stat_compare_means(comparisons = comp)

m_sub <- subset(mdf, module == "M1" | module == "M2" | module == "M3" | module == "M4")

ggboxplot(m_sub, x="Treatment", y="RelativeAbundance", fill = "Phylum", facet.by = "module", scales = "free_y")

bz <- subset(m_sub, Treatment == "bzwf" | Treatment == "bznf")
sarc <- subset(m_sub, Treatment == "sarcwf" | Treatment == "sarcnf")

ggboxplot(bz, x="Treatment", y="RelativeAbundance", fill = "Phylum", facet.by = "module")
ggboxplot(sarc, x="Treatment", y="RelativeAbundance", fill = "Phylum", facet.by = "module")

ggbarplot(m_sub, x="Treatment", y="RelativeAbundance", fill = "Phylum", facet.by = "module", scales = "free_y")

# Make a microeco object with individual Modules
m1_asv <- as.data.frame(unique(m1$name))
m2_asv <- as.data.frame(unique(m2$name))
m3_asv <- as.data.frame(unique(m3$name))
m4_asv <- as.data.frame(unique(m4$name))

rownames(m1_asv) <- m1_asv[,1]
rownames(m2_asv) <- m2_asv[,1]
rownames(m3_asv) <- m3_asv[,1]
rownames(m4_asv) <- m4_asv[,1]

m1_meco <- clone(m16)
m2_meco <- clone(m16)
m3_meco <- clone(m16)
m4_meco <- clone(m16)

m1_meco$otu_table %<>% subset(rownames(m1_meco$otu_table) %in% rownames(m1_asv))
m2_meco$otu_table %<>% subset(rownames(m2_meco$otu_table) %in% rownames(m2_asv))
m3_meco$otu_table %<>% subset(rownames(m3_meco$otu_table) %in% rownames(m3_asv))
m4_meco$otu_table %<>% subset(rownames(m4_meco$otu_table) %in% rownames(m4_asv))

m1_meco$tidy_dataset()
m2_meco$tidy_dataset()
m3_meco$tidy_dataset()
m4_meco$tidy_dataset()

m1_meco$cal_betadiv()
m2_meco$cal_betadiv()
m3_meco$cal_betadiv()
m4_meco$cal_betadiv()

m1_meco$sample_table
m1_meco$sample_table$Loc_Fus <- paste0(m1_meco$sample_table$Location, "_", m1_meco$sample_table$FusariumPresence)

tm1 <- trans_abund$new(dataset = m1_meco, taxrank = "Phylum", groupmean = "Loc_Fus")
tm1$plot_bar()

bz_m1 <- clone(m1_meco)
sarc_m1 <- clone(m1_meco)

bz_m1$sample_table %<>% subset(Location == "Bozeman")
sarc_m1$sample_table %<>% subset(Location == "Huntley")

bz_m1$tidy_dataset()
sarc_m1$tidy_dataset()

bz_tm1 <- trans_abund$new(dataset = bz_m1, taxrank = "Phylum", ntaxa = 10)
bz_tm1$plot_box(group = "FusariumPresence", xtext_angle = 30) 

dif_bz_m1 <- trans_diff$new(dataset = bz_m1, method = "lefse", group = "Fus_cat", taxa_level = "Genus")
dif_bz_m1$plot_diff_abund()

dif_sarc_m1 <- trans_diff$new(dataset = sarc_m1, method = "lefse", group = "Fus_cat", taxa_level = "Genus")
dif_sarc_m1$plot_diff_abund()

# show top 10 taxa at Family level
tfm1 <- trans_abund$new(dataset = m1_meco, taxrank = "Family", ntaxa = 10)
tfm1$plot_box(group = "Loc_Fus", xtext_angle = 30)

dif_m1 <- trans_diff$new(dataset = m1_meco, method = "lefse", group = "Loc_Fus", taxa_level = "Phylum")
# see t1$res_diff for the result
# From v0.8.0, threshold is used for the LDA score selection.
dif_m1$plot_diff_bar(threshold = 4)
dif_m1$plot_diff_abund(use_number = 1:10)

# Correlation
m1_sd <- as.data.frame(m1_meco$sample_table)
m2_sd <- as.data.frame(m2_meco$sample_table)
m3_sd <- as.data.frame(m3_meco$sample_table)
m4_sd <- as.data.frame(m4_meco$sample_table)

m1_otu <- as.data.frame(m1_meco$otu_table)
m2_otu <- as.data.frame(m2_meco$otu_table)
m3_otu <- as.data.frame(m3_meco$otu_table)
m4_otu <- as.data.frame(m4_meco$otu_table)

m1_otu <- t(m1_otu)
m2_otu <- t(m2_otu)
m3_otu <- t(m3_otu)
m4_otu <- t(m4_otu)

m1_otu <- as.data.frame(m1_otu)
m2_otu <- as.data.frame(m2_otu)
m3_otu <- as.data.frame(m3_otu)
m4_otu <- as.data.frame(m4_otu)

m1_otu$Sum <- rowSums(m1_otu[,c(1:188)])
m2_otu$Sum <- rowSums(m2_otu[,c(1:156)])
m3_otu$Sum <- rowSums(m3_otu[,c(1:129)])
m4_otu$Sum <- rowSums(m4_otu[,c(1:109)])

m1_otu <- m1_otu[,c(189,1:188)]
m2_otu <- m2_otu[,c(157,1:156)]
m3_otu <- m3_otu[,c(130,1:129)]
m4_otu <- m4_otu[,c(110,1:109)]

m1_otu <- m1_otu[,c(1,2)]
m2_otu <- m2_otu[,c(1,2)]
m3_otu <- m3_otu[,c(1,2)]
m4_otu <- m4_otu[,c(1,2)]

m1_df <- cbind(m1_sd, m1_otu)
m2_df <- cbind(m2_sd, m2_otu)
m3_df <- cbind(m3_sd, m3_otu)
m4_df <- cbind(m4_sd, m4_otu)

m1_df <- m1_df[,c(17,35)]
m2_df <- m2_df[,c(17,35)]
m3_df <- m3_df[,c(17,35)]
m4_df <- m4_df[,c(17,35)]

ggscatter(m1_df, x = "Fusarium", y = "Sum", add = "reg.line") + ylab("ASV Abundance") + xlab("Number of Fusarium Isolates") +
  stat_cor(label.x = 3, label.y = 130000)

ggscatter(m2_df, x = "Fusarium", y = "Sum", add = "reg.line") + ylab("ASV Abundance") + xlab("Number of Fusarium Isolates") +
  stat_cor(label.x = 3, label.y = 30000)

ggscatter(m3_df, x = "Fusarium", y = "Sum", add = "reg.line") + ylab("ASV Abundance") + xlab("Number of Fusarium Isolates") +
  stat_cor(label.x = 3, label.y = 75000)

ggscatter(m4_df, x = "Fusarium", y = "Sum", add = "reg.line") + ylab("ASV Abundance") + xlab("Number of Fusarium Isolates") +
  stat_cor(label.x = 3, label.y = 25000)

  
# Hub taxa calculated with Cytohubba in Cytoscape and .csvs imported into R for analysis
bzb_hub <- read.csv("Bozeman Fusarium_MCC_top20 default node.csv")
bznb_hub <- read.csv("Bozeman No Fusarium_MCC_top20 default node.csv")
sarcb_hub <- read.csv("Huntley Fusarium_MCC_top20 default node.csv")
sarcnb_hub <- read.csv("Huntley No Fusarium_MCC_top20 default node.csv")

bzb_hub$Location <- "Bozeman"
bznb_hub$Location <- "Bozeman"
sarcb_hub$Location <- "Huntley"
sarcnb_hub$Location <- "Huntley"

bzb_hub$Fusarium <- "Fusarium"
bznb_hub$Fusarium <- "No Fusarium"
sarcb_hub$Fusarium <- "Fusarium"
sarcnb_hub$Fusarium <- "No Fusarium"

hub_df <- rbind(bzb_hub, bznb_hub, sarcb_hub, sarcnb_hub)
write.csv(hub_df, "Hub Taxa Loc Fusarium.csv")

unique(hub_df$Phylum)
hub_df$Phylum <- factor(hub_df$Phylum, levels = c("Acidobacteria", "Actinobacteria",  "Bacteroidetes", "Firmicutes", "Armatimonadetes", "Proteobacteria", "Cyanobacteria"))
ggbarplot(hub_df, x = "Fusarium", y= "RelativeAbundance", fill = "Phylum", facet.by = "Location") + ylab("Relative Abundance (%)")

unique(hub_df$Family)
ggbarplot(hub_df, x = "Fusarium", y= "RelativeAbundance", fill = "Family", facet.by = "Location") + ylab("Relative Abundance (%)")

# Venn Diagram
library(ggvenn)

# Phylum
x2 <- list('Bozeman Fusarium' = bzb_hub$Phylum, 'Bozeman No Fusarium' = bznb_hub$Phylum, "Huntley No Fusarium" = sarcnb_hub$Phylum, 'Huntley Fusarium' = sarcb_hub$Phylum)

pdf("Figures/Venn Phylum Hub.pdf", width = 10, height = 8)
ggvenn(
  x2, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#CD80C2FF" ,"#ED73C2FF"),
  stroke_size = 0, set_name_size = 4, show_outside = "always"
)
dev.off() 

# Family
x3 <- list('Bozeman Fusarium' = bzb_hub$Family, 'Bozeman No Fusarium' = bznb_hub$Family, "Huntley No Fusarium" = sarcnb_hub$Family, 'Huntley Fusarium' = sarcb_hub$Family)

pdf("Figures/Venn Family Hub.pdf", width = 10, height = 8)
ggvenn(
  x3, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#CD80C2FF" ,"#ED73C2FF"),
  stroke_size = 0, set_name_size = 4, show_outside = "always"
)
dev.off() 

# Genus
x4 <- list('Bozeman Fusarium' = bzb_hub$Genus, 'Bozeman No Fusarium' = bznb_hub$Genus, "Huntley No Fusarium" = sarcnb_hub$Genus, 'Huntley Fusarium' = sarcb_hub$Genus)

pdf("Figures/Venn Genus Hub.pdf", width = 10, height = 8)
ggvenn(
  x4, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#CD80C2FF" ,"#ED73C2FF"),
  stroke_size = 0, set_name_size = 4, show_outside = "always"
)
dev.off() 

# Differential abundance analysis of hub taxa
library(ggnested)
m16 <- readRDS("../m16S_combined.RDS")
dim(m16$tax_table)

m16S_sub <- clone(m16)

m16S_sub$tax_table <- m16S_sub$tax_table[rownames(m16S_sub$tax_table) %in% hub_df$id, ]
dim(m16$tax_table)
dim(m16S_sub$tax_table)
m16S_sub$tidy_dataset()


t1 <- trans_abund$new(dataset = m16S_sub, taxrank = "Genus", ntaxa = 76, high_level = "Phylum")
t1$plot_box(group = "Loc_Yr", xtext_angle = 30)
t1$plot_bar(facet = c("Loc_Yr", "FusariumPresence"), ggnested = TRUE)
t1$plot_heatmap(facet = "Loc_Yr", xtext_keep = FALSE, withmargin = FALSE, plot_breaks = c(0.01, 0.1, 1, 10))

m16S_sub$sample_table
t2 <- trans_diff$new(dataset = m16S_sub, method = "anova", group = "Loc_Yr", taxa_level = "Genus")
t2$plot_diff_abund(plot_type = "barerrorbar", errorbar_addpoint = FALSE)

