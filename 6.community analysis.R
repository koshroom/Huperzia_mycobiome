setwd("/Users/ko/Work/AS_project/Huperzia")
library("ape")
library("vegan")
library("tidyr")
#library("metacoder")
library(phyloseq)
library("ggplot2")
#library("indicspecies")
library(EcolUtils)
library(ranacapa)
library(DESeq2)
library("knitr")
library("BiocStyle")
library("dada2")
library("DECIPHER")
library(microbiome)
library(RColorBrewer)
library(dplyr)
library(purrr)
library("phangorn")
library("Biostrings")
library("lulu")
library("metagMisc")
library("MiscMetabar")
library(FUNGuildR)
library(tidyverse)
library(moments)
library(ggpubr)
library(microbiomeMarker)
library(multcompView)

## Load taxa produced by dada2##
taxa <- readRDS(file = "taxa_euc.rds")
tax_R <- gsub(".__", "", taxa) ## If having problem with tax rank, skip this

meta <-read.delim("hup_metadata.txt", row.names = 1)

## Load seqtab.nochim produced by dada2##
seqtab.nochim <- readRDS(file = "seqtab.nochim.rds")
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
row.names(seqtab.nochim) <- gsub("LGL21_EO01_PCR_R_(bc[0-9]+)_.+","\\1", row.names(seqtab.nochim))

#########################
##### Start Phyloseq#####

# ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
#                sample_data(meta), 
#                tax_table(tax_R), phy_tree(fit$tree))
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(meta),
               tax_table(tax_R))
ps
#unique(data.frame(tax_table(ps))$Kingdom) #Just double check if some ASV belong to plants

dna <- Biostrings::DNAStringSet(taxa_names(ps)) ## Add refsequences to phyloseq object
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps


luall <- lulu_phyloseq(
  ps,
  nproc = 1,
  id = 0.84,
  vsearchpath = "/./Users/ko/software/vsearch-2.21.1-macos-x86_64/bin/vsearch",
  verbose = FALSE,
  clean_physeq = FALSE
)


luall_ps <- luall$new_physeq
ok_seq <- names(refseq(luall_ps))[(width(refseq(luall_ps))) > 150]  ## Filter short ASVs
luall_ps <- prune_taxa(ok_seq, luall_ps)
#saveRDS(luall_ps, "luall_ps.RDS")

## To measure plant:fungal ratio ####
luall_ps_h <-subset_samples(luall_ps, species == "Huperzia_selago")
luall_psf_h <- subset_taxa(luall_ps_h, Kingdom == "Fungi")

fungalkingdom <- otu_table(tax_glom(luall_psf_h, taxrank = "Kingdom"))
phylo_lyco <- subset_taxa(luall_ps_h, Phylum == "Lycopodiophyta")
lycophylum <- otu_table(tax_glom(phylo_lyco, taxrank = "Phylum"))
read_ratio <- merge(fungalkingdom, lycophylum, by = "row.names")
colnames(read_ratio) <- c("sample", "fungi", "plant")
row.names(read_ratio) <- read_ratio$sample
read_ratio[read_ratio$plant == 0, "plant"] <- 1
read_ratio$ratio <- with(read_ratio, fungi/plant)
concate <- merge(read_ratio, as.data.frame(sample_data(luall_psf_h)), by = "row.names")
concate[concate$tissue == "leaf", "tissue"] <- "standard_leaf"
plot(concate)

png("fungi_plant_ratio.png", width = 3400, height = 3000, res = 900)
ggplot(concate, aes(x = factor(tissue,  level=c("bulbil", "shoot", "young_leaf","standard_leaf", "stem","sporangium")), y = ratio)) +
  geom_point() +  xlab("tissue") +
  ylab("fungi/plant ratio") +
  geom_boxplot() +
  theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off()

####
psf <- subset_taxa(ps, Kingdom == "Fungi") ## From here include fungi only
ok_seq <- names(refseq(psf))[(width(refseq(psf))) > 150]  ## Filter short ASVs
psf <- prune_taxa(ok_seq, psf)
psf

### Lulu assignment ####
lu <- lulu_phyloseq(
  psf,
  nproc = 1,
  id = 0.84,
  vsearchpath = "/./Users/ko/software/vsearch-2.21.1-macos-x86_64/bin/vsearch",
  verbose = FALSE,
  clean_physeq = FALSE
)

lu_psf <- lu$new_physeq
ok_seq <- names(refseq(lu_psf))[(width(refseq(lu_psf))) > 150]  ## Filter short ASVs
lu_psf <- prune_taxa(ok_seq, lu_psf)

writeXStringSet(DNAStringSet(refseq(lu_psf)), "lu_refseq.fasta", append=FALSE, #Prep for vseach for Lulu #/Users/ko/software/vsearch-2.21.1-macos-x86_64/bin
                compress=FALSE, compression_level=NA, format="fasta") #Will be used to cluster with HupA producing fungi


#### Funguild assignment ####
guild <- as_tibble(as.data.frame(tax_table(psf)), rownames = NA) %>%  
  rownames_to_column %>% ## Need to figure out how to keep rownames
  unite(Taxonomy, "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", sep = ";", remove = FALSE) %>% 
  funguild_assign

guild_test <- as.matrix(guild)
row.names(guild_test) <- guild_test[,1]
guild_test <- guild_test[,-1]
dim(tax_table(lu_psf))
head(tax_table(guild_test))
tax_table(psf) <- tax_table(guild_test)


guild <- as_tibble(as.data.frame(tax_table(lu_psf)), rownames = NA) %>%  
  rownames_to_column %>% ## Need to figure out how to keep rownames
  unite(Taxonomy, "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", sep = ";", remove = FALSE) %>% 
  funguild_assign

guild_test <- as.matrix(guild)
row.names(guild_test) <- guild_test[,1]
guild_test <- guild_test[,-1]
dim(tax_table(lu_psf))
head(tax_table(guild_test))
tax_table(lu_psf) <- tax_table(guild_test)


##### Color define ####
################
phy_bar_color <- c(Ascomycota = "red", Basidiomycota = "blue", Mucoromycota = "#66CC99", `NA` = "white")

tissuecol_all <- c(bulbil = "#C77CFF", shoot = "#87CEEB", young_leaf = "#66CC99", standard_leaf = "seagreen4", 
                   stem = "#FF7F50", sporangium ="#DAA520")

# tissuecol_all <- c(standard_leaf="seagreen4", shoot = "#87CEEB",
#                    stem = "#FF7F50", bulbil = "#C77CFF",
#                    sporangium ="#DAA520", young_leaf = "#66CC99", root_cluster = "#663300", root = "#330000")

tissuecol_h <- c(standard_leaf = "seagreen4", shoot = "#87CEEB",
                 stem = "#FF7F50", bulbil = "#C77CFF",
                 sporangium ="#DAA520", young_leaf = "#66CC99")

tissue_stem_leaf <- c(standard_leaf="seagreen4", stem = "#FF7F50")
tissue_stem_leaf1 <- c(standard_leaf="seagreen4", stem = "#FF7F50", young_leaf = "#66CC99")

phylum_col <- c(Ascomycota = "plum", Basidiomycota = "gold3", Mucoromycota = "hotpink4")
##

# levels(sample_data(psf_p)$tissue) = c("bulbil", "shoot", "young_leaf",
#                                        "standard_leaf", "stem","sporangium") #do we need this after adding the "breaks comment"?
#                                        
# 



# psf <- readRDS("lu_psf.rds") #important
sample_data(psf)$tissue[sample_data(psf)$tissue == "leaf"] <- "standard_leaf" #change the metadata of original "leaf" to "old_leaf". added 2023
psf <- subset_samples(psf, !tissue %in% c("root", "root_cluster")) #discard root related samples from Lycopodium
psf <- subset_taxa(psf, !is.na(Phylum)) ## Remove phylum = NA
## proportion ##
filter1 <- phyloseq::genefilter_sample(psf, filterfun_sample(function(x) x >= 1),
                                      A=1)
psf_filtered1 <- prune_taxa(filter1, psf) #just to check the number of taxa remained
psf <- psf_filtered1 


filter <- phyloseq::genefilter_sample(psf, filterfun_sample(function(x) x >= 2),
                                      A=2)
psf_filtered <- prune_taxa(filter, psf)
psf_filtered 
psf_p <- transform_sample_counts(psf_filtered, function(x){x / sum(x)}) 
psf_p


#### Create Rarified table with singleton###
nreads = sort(sample_sums(psf)) # if need the lowest read number for a sample to determine the rarification number
set.seed(1111)
mean(nreads)
sd(nreads)

psf_R = rarefy_even_depth(psf, sample.size = nreads[[1]])

nread = sort(sample_sums(psf)) # if need the lowest read number for a sample to determine the rarification number
nreads[[1]]


title = "Sum of reads for each sample, physeq_R"
plot(sort(sample_sums(psf_R), TRUE), type = "h", main = title, ylab = "reads", 
     ylim = c(0, 20000)) ## to check if the rarefaction did work



## subset data ## 
###

psf_h <- subset_samples(psf, species == "Huperzia_selago")
psf_l <- subset_samples(psf, species == "Lycopodium_complanatum")
filter_h <- phyloseq::genefilter_sample(psf_h, filterfun_sample(function(x) x >= 2),
                                                                              A=2)

psf_h_youngleaf <- subset_samples(psf_h, tissue == "young_leaf")
#sample_data(psf_h)$tissue[sample_data(psf_h)$tissue == "leaf"] <- "standard_leaf" 
#psf_l_leaf <- subset_samples(psf_l, tissue == "standard_leaf")

#psf_l_h_younleaf <- merge_phyloseq(psf_h_youngleaf, psf_l_leaf)

psf_h_filtered <- prune_taxa(filter_h, psf_h)
# psf_h_oleaf_filtered <- subset_samples(psf_h_filtered, tissue == "old_leaf")
# psf_h_yleaf_filtered <- subset_samples(psf_h_filtered, tissue == "young_leaf")
# psf_h_bothleaf_filtered <- merge_phyloseq(psf_h_oleaf_filtered, psf_h_yleaf_filtered)


psf_leaf <- subset_samples(psf, tissue %in% c("leaf", "standard_leaf", "young_leaf"))
psf_stem <- subset_samples(psf, tissue == "stem")
psf_stleaf <- subset_samples(psf, tissue == "standard_leaf")


psf_p_h <- subset_samples(psf_p, species == "Huperzia_selago")
psf_p_l <- subset_samples(psf_p, species == "Lycopodium_complanatum")
psf_R_h <- subset_samples(psf_R, species == "Huperzia_selago")
psf_R_leaf <- subset_samples(psf_R, tissue %in% c("leaf", "standard_leaf", "young_leaf"))
psf_R_stem <- subset_samples(psf_R, tissue == "stem")
psf_p_leaf <- subset_samples(psf_p, tissue %in% c("leaf", "standard_leaf", "young_leaf"))
psf_p_stem <- subset_samples(psf_p, tissue == "stem")
psf_p_leafstem <- merge_phyloseq(psf_p_leaf, psf_p_stem)
psf_R_leafstem <- merge_phyloseq(psf_R_leaf, psf_R_stem)

## Figures ####

psf_p.NMDS = ordinate(psf_p, "NMDS", "bray", trymax = 1000) #
ptitle = "brayNMDS taxa >2 present, stress 0.11"

png("NMDS_bray_all.png", width=4400, height=3000,res=300)
p = plot_ordination(psf_p, psf_p.NMDS, shape = "species", 
                    color = "tissue") +stat_ellipse()
p=p + geom_point(size = 9)+ 
  ggtitle(ptitle)+ 
  scale_color_manual(values=tissuecol_all, breaks = c("bulbil", "shoot", "young_leaf", "standard_leaf", "stem","sporangium"))+ 
                             theme(text= element_text(size=30),
                              axis.text = element_text(size=40), 
                            panel.background=element_rect(fill="white",colour="black"))
p
dev.off()


#
psf_p_h.NMDS = ordinate(psf_p_h, "NMDS", "bray", trymax = 10000) #
ptitle = "brayNMDS taxa >2 present"

#png("NMDS_h_all.png", width=4400, height=3000,res=300)
pdf("NMDS_h_all.pdf", width=10, height=8)
p = plot_ordination(psf_p_h, psf_p_h.NMDS, 
                    color = "tissue") + stat_ellipse()
p=p + geom_point(size = 7)+ ggtitle(ptitle)+ scale_color_manual(values=tissuecol_h, breaks = c("bulbil", "shoot", "young_leaf", "standard_leaf", "stem","sporangium")) + theme(text= element_text(size=30), 
                                                   axis.text = element_text(size=40), 
                                                   panel.background=element_rect(fill="white",colour="black"))
p
dev.off()


psf_p_leafstem.NMDS = ordinate(psf_p_leafstem, "NMDS", "bray", trymax = 1000) #
ptitle = "brayNMDS taxa >2, proportion 0.10"

#png("NMDS_p_leafstem_all.png", width=4400, height=3000,res=300)
pdf("NMDS_p_leafstem_all.pdf", width=12, height=8)
p = plot_ordination(psf_p_leafstem, psf_p_leafstem.NMDS, 
                    color = "tissue", shape = "species") + stat_ellipse()
p=p + geom_point(size = 9)+ ggtitle(ptitle)+ scale_color_manual(values=tissue_stem_leaf1, breaks = c("young_leaf", "standard_leaf", "stem")) + theme(text= element_text(size=30), 
                                                                                            axis.text = element_text(size=40), 
                                                                                            panel.background=element_rect(fill="white",colour="black"))
p
dev.off()





## Richness ##
png("psf_R_h_rich.png", width=3200, height=2400, res=250)
plot_richness(psf_R_h, x = "tissue") + geom_boxplot()
dev.off()


rich_psf_R_h_tissue <- data.frame(
  "Observed" = phyloseq::estimate_richness(psf_R_h, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(psf_R_h, measures = "Shannon"),
  "tissue" = phyloseq::sample_data(psf_R_h)$tissue)

shapiro.test(rich_psf_R_h_tissue$Observed)
shapiro.test(rich_psf_R_h_tissue$Shannon)


moments::skewness(rich_psf_R_h_tissue$Observed, na.rm = TRUE) #https://bookdown.org/wadetroberts/r-you-ready-for-r/anova.html#packages-needed-for-anova
moments::kurtosis(rich_psf_R_h_tissue$Observed, na.rm = TRUE)
moments::skewness(rich_psf_R_h_tissue$Shannon, na.rm = TRUE) #https://bookdown.org/wadetroberts/r-you-ready-for-r/anova.html#packages-needed-for-anova
moments::kurtosis(rich_psf_R_h_tissue$Shannon, na.rm = TRUE)
# ggpubr::ggdensity(mydata$intvar, fill = "lightgray")
# ggpubr::ggqqplot(mydata$intvar)
kruskal.test(Observed ~ tissue, data = rich_psf_R_h_tissue)
kruskal.test(Shannon ~ tissue, data = rich_psf_R_h_tissue)




#png("Observed number of species of Huperzia selago tissues.png", width=3000, height=2500, res=600)
pdf("Observed number of species of Huperzia selago tissues.pdf", width= 7, height= 6)
p <- rich_psf_R_h_tissue %>%
  ggplot(aes(x=reorder(tissue, Observed, na.rm = TRUE), y = Observed)) +
  geom_boxplot() +
  geom_point()+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1))
  
p+ labs(y="Observed number of species", x="tissue", 
        subtitle="Observed number of species of Huperzia selago tissues")
dev.off()


#png("Shannon diversity of Huperzia selago tissues.png", width=3000, height=2500, res=600)
pdf("Shannon diversity of Huperzia selago tissues.pdf", width=7, height=6)
p1 <- rich_psf_R_h_tissue %>%
  ggplot(aes(x=reorder(tissue, Shannon, na.rm = TRUE), y = Shannon)) +
  geom_boxplot() +
  geom_point()+
  theme_bw()+ theme(axis.text.x=element_text(angle=45, hjust=1))
  
p1 +  labs(y="Shannon diversity", x="tissue", 
       subtitle="Shannon diversity of Huperzia selago tissues")
dev.off()

# png("test")
# par(mfrow = c(2,2))
# plot(p)
# plot(p1)
# dev.off()

####### TRy to add multicom letters##### I added them manually # Observed is normally distributed while Shannon was not (not after log transformed either)

dunn.test::dunn.test(rich_psf_R_h_tissue$Observed, rich_psf_R_h_tissue$tissue, method = "holm")
test1 <- kruskal.test(Observed ~ tissue, data = rich_psf_R_h_tissue)
test2 <-dunnTest(Observed ~ tissue, data = rich_psf_R_h_tissue)
multcompLetters4(test1, test2)#notworking
anova <- aov(Observed ~ tissue, data = log_rich_psf_R_h_tissue)
Tukey <- TukeyHSD(anova)
cld <- multcompLetters4(anova, Tukey)
# # Table with the mean, the standard deviation and the letters indications significant differences for each treatment
# dt <- group_by(rich_psf_R_h_tissue, tissue) %>%
#   summarise(o_mean=mean(Observed), sd=sd(Observed)) %>%
#   arrange(desc(o_mean))
# cld <- as.data.frame.list(cld$`factor(tissue)`)
# dt$Tukey <- cld$Letters
# 
# print(dt)
# # barplot
# ggplot(dt, aes(x = factor(conc), y = uptake_mean, fill = Treatment)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   geom_errorbar(aes(ymax = uptake_mean + sd, ymin = uptake_mean - sd),
#                 position = position_dodge(0.9), width = 0.25, color = "Gray25") +
#   xlab(expression(CO[2]~Concentration~'('~mL~L^-1~')')) +
#   ylab(expression(CO[2]~Uptake~'('~µmol~m^2~s^-1~')')) +
#   theme_few() +
#   theme(legend.position = c(0.58, 0.8)) +
#   scale_fill_brewer(palette = "Greens") +
#   facet_grid(.~Type, labeller = label_both) +
#   geom_text(aes(label=Tukey, y = uptake_mean + sd + 2), size = 3, color = "Gray25",
#             show.legend = FALSE,
#             position = position_dodge(0.9)) +
#   ylim(0,50)









rich_psf_R_leafstem <- data.frame(
   "Observed" = phyloseq::estimate_richness(psf_R_leafstem, measures = "Observed"),
   "Shannon" = phyloseq::estimate_richness(psf_R_leafstem, measures = "Shannon"),
   "tissue" = phyloseq::sample_data(psf_R_leafstem)$tissue,
   "species" = phyloseq::sample_data(psf_R_leafstem)$species)

png("psf_R_leafstem_rich.png", width=3200, height=2400, res=650) #likely wrong, probably combine shannon & observed
plot_richness(psf_R_leafstem, x = "tissue") + geom_boxplot() + facet_wrap(~species)+
  theme_bw()+ theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off() 

# png("Observed number of species of leaf_stem tissues.png", width=2800, height=2500, res=600)
#  rich_psf_R_h_tissue %>%
#    ggplot(aes(x=reorder(tissue, Observed, na.rm = TRUE), y = Observed)) +
#    geom_boxplot() +
#    geom_point()+
#    theme_bw()+
#    labs(y="Observed number of species", x="tissue", 
#         subtitle="Observed number of species of Huperzia selago tissues") +
#    facet_wrap(~species)
#  dev.off()
# 

###################
## Stack barplot ##
################### # https://micca.readthedocs.io/en/latest/phyloseq.html

psf_h_class = tax_glom(psf_h, taxrank="Class", NArm=FALSE)
psf_p_h_phylum = tax_glom(psf_p_h, taxrank="Phylum", NArm=FALSE)
psf_p_h_class = tax_glom(psf_p_h, taxrank="Class", NArm=FALSE)


psf_p_h <- readRDS(psf_p_h)

psf_p_h_phylum  <- psf_p_h %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Phylum)   

png("barplot_h_phy.png", width=8000, height=3000,res=800)
ggplot(psf_p_h_phylum, aes(x= Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") +
  ylab("Relative Abundance") +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  scale_color_manual(values =phy_bar_color)
dev.off()




psf_p_h_class  <- psf_p_h %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at Family level
  # transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Class)   

psf_p_l_class  <- psf_p_l %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at Family level
  # transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Class)   


psf_p_class  <- psf_p %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at Family level
  # transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Class)   



png("barplot_huperzia_non100.png", width=5000, height=3000,res=350)
title = "Huperzia_class"
plot_bar(psf_h_class, fill="Class") + facet_wrap(~tissue, scales= "free_x", nrow=1)
dev.off()


png("barplot_huperzia_p.png", width=3000, height=3000,res=350)
title = "Huperzia_phylum"
plot_bar(psf_p_h_phylum, fill="Phylum") + facet_wrap(~tissue, scales= "free_x", nrow=1)
dev.off()


pdf("barplot_h_class.pdf", width=10, height=6)
#png("barplot_h_class.png", width=8000, height=3000,res=800)
ggplot(psf_p_h_class, aes(x= Sample, y = Abundance, fill = Class)) + 
  geom_bar(stat = "identity") +
  #  scale_fill_manual(values=taxcol)+
  ylab("Relative Abundance (Class > 0.1%) \n") +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  #ggtitle("Class level ITS1") +
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~ factor(tissue, levels = c("bulbil", "shoot", "young_leaf", "standard_leaf", "stem","sporangium")), scales= "free_x", nrow=1)
dev.off()


png("barplot_l_class.png", width=4000, height=3000,res=800)
ggplot(psf_p_l_class, aes(x= Sample, y = Abundance, fill = Class)) + 
  geom_bar(stat = "identity") +
  #  scale_fill_manual(values=taxcol)+
  ylab("Relative Abundance (Class > 0.1%) \n") +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  #ggtitle("Class level ITS1") +
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~ factor(tissue, levels = c("standard_leaf", "stem")), scales= "free_x", nrow=1)
dev.off()


png("barplot_all_class.png", width=12000, height=3000,res=800)
ggplot(psf_p_class, aes(x= Sample, y = Abundance, fill = Class)) + 
  geom_bar(stat = "identity") +
  #  scale_fill_manual(values=taxcol)+
  ylab("Relative Abundance (Class > 0.1%) \n") +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  #ggtitle("Class level ITS1") +
  scale_fill_brewer(palette = "Paired")+
  facet_wrap(species ~ factor(tissue, levels = c("bulbil", "shoot", "young_leaf", "standard_leaf", "stem","sporangium")), scales= "free_x", nrow=1)
dev.off()


################
#### DEseq2 ####
################ https://joey711.github.io/phyloseq-extensions/DESeq2.html

##### read HupA ASV #####
vsearch <- read.delim("/Users/ko/Work/AS_project/Huperzia/share/nov2022/vsearch_ASV_97_whuperzia_culture.txt", header = F)
dim(vsearch)
geneious <- read.delim("/Users/ko/Work/AS_project/Huperzia/share/nov2022/Geneious_230\ documents\ from\ Huperzia.tsv", header = F)
dim(geneious)
hupAF <- unique(c(vsearch[[1]], geneious[[1]]))
write.table(hupAF, "hupA_producing_ASV.txt", sep = "\t", quote = F, row.names = F)

tax_table(psf)
hupAF <- data.frame(hupAF)
row.names(hupAF) <- hupAF$hupAF
hupAF <- merge(hupAF, data.frame(tax_table(psf)), by = "row.names", all.x =TRUE)
hupAF <- hupAF[,-1]
write.table(hupAF, "hupA_producing_ASV_wtaxainfo.txt", sep = "\t", quote = F, row.names = F)
## Compare the leaf mycobiome of two species ##
###############################################

psf_leaf
psf_leaf_de = phyloseq_to_deseq2(psf_leaf, ~ species)
psf_leaf_de = DESeq(psf_leaf_de, test="Wald", fitType="parametric") 
res = results(psf_leaf_de, cooksCutoff = FALSE)
summary(res)
alpha = 0.01

#saveRDS(res, "res.RDS")
#res <- readRDS("res_LleafvsH_oldleaf.RDS")

# upregulated in Huperzia
sigtab_Hup = res[which(res$padj < alpha & res$log2FoldChange<0 & abs(res$log2FoldChange)>2), ]
sigtab_Hup = cbind(as(sigtab_Hup, "data.frame"), as(tax_table(psf_leaf)[rownames(sigtab_Hup), ], "matrix"))
head(sigtab_Hup)
dim(sigtab_Hup)

sigtab_Lup = res[which(res$padj < alpha & res$log2FoldChange>0 & abs(res$log2FoldChange)>2), ]
sigtab_Lup = cbind(as(sigtab_Lup, "data.frame"), as(tax_table(psf_leaf)[rownames(sigtab_Lup), ], "matrix"))
head(sigtab_Lup)
dim(sigtab_Lup)

sigtab_Hup_hupAF <- sigtab_Hup[row.names(sigtab_Hup) %in% hupAF,]
dim(sigtab_Hup_hupAF)
sigtab_Lup_hupAF <- sigtab_Lup[row.names(sigtab_Lup) %in% hupAF,]
dim(sigtab_Lup_hupAF)

cor_sigtab_Hup_hupAF <- sigtab_Hup[row.names(sigtab_Hup) %in% names(ASV_hupA_0.4),]
dim(cor_sigtab_Hup_hupAF)
cor_sigtab_Lup_hupAF <- sigtab_Lup[row.names(sigtab_Lup) %in% names(ASV_hupA_0.4),]
dim(cor_sigtab_Lup_hupAF)

## upregulated in H
x = tapply(sigtab_Hup$log2FoldChange, sigtab_Hup$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_Hup$Phylum = factor(as.character(sigtab_Hup$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab_Hup$log2FoldChange, sigtab_Hup$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab_Hup$Genus = factor(as.character(sigtab_Hup$Genus), levels=names(x))
# 
# theme_set(theme_bw()) #delete?
# scale_fill_discrete <- function(palname = "Set1", ...) {
#   scale_fill_brewer(palette = palname, ...)
# }
#sigtab_Hup <- readRDS("sigtab_Hup.RDS")
#png("logfold_leaf_Hup_genus031923.png", width=5000, height=3000,res=550)
pdf("logfold_leaf_Hup_genus031923.pdf", width=10, height=6)
ggplot(sigtab_Hup, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=5) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  scale_color_manual(values = phylum_col)
dev.off()


#up regulated in L
# Phylum order
x = tapply(sigtab_Lup$log2FoldChange, sigtab_Lup$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_Lup$Phylum = factor(as.character(sigtab_Lup$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab_Lup$log2FoldChange, sigtab_Lup$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab_Lup$Genus = factor(as.character(sigtab_Lup$Genus), levels=names(x))


# theme_set(theme_bw()) #delete?
# scale_fill_discrete <- function(palname = "Set1", ...) {
#   scale_fill_brewer(palette = palname, ...)
# }

#sigtab_Lup <- readRDS("sigtab_Lup.RDS")
#png("logfold_leaf_Lup_genus_031923.png", width=5000, height=3000,res=550)
pdf("logfold_leaf_Lup_genus_031923.pdf", width=10, height=6)
ggplot(sigtab_Lup, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=5) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  scale_color_manual(values = phylum_col) 
dev.off()

saveRDS(res, "res_LleafvsH_oldleaf.RDS")
saveRDS(sigtab_Hup, "sigtab_Hup.RDS")
saveRDS(sigtab_Lup, "sigtab_Lup.RDS")
saveRDS(sigtab_Hup_hupAF, "sigtab_Hup_hupAF.RDS")
saveRDS(sigtab_Lup_hupAF, "sigtab_Lup_hupAF.RDS")



##### 
psf_stleaf ## standard leaf only

psf_stleaf
psf_stleaf_de = phyloseq_to_deseq2(psf_stleaf, ~ species)
psf_stleaf_de = DESeq(psf_stleaf_de, test="Wald", fitType="parametric") 
res = results(psf_stleaf_de, cooksCutoff = FALSE)
summary(res)
alpha = 0.01

# upregulated in Huperzia
sigtab_Hup = res[which(res$padj < alpha & res$log2FoldChange<0 & abs(res$log2FoldChange)>2), ]
sigtab_Hup = cbind(as(sigtab_Hup, "data.frame"), as(tax_table(psf_leaf)[rownames(sigtab_Hup), ], "matrix"))
head(sigtab_Hup)
dim(sigtab_Hup)

sigtab_Lup = res[which(res$padj < alpha & res$log2FoldChange>0 & abs(res$log2FoldChange)>2), ]
sigtab_Lup = cbind(as(sigtab_Lup, "data.frame"), as(tax_table(psf_leaf)[rownames(sigtab_Lup), ], "matrix"))
head(sigtab_Lup)
dim(sigtab_Lup)


sigtab_Hup_hupAF <- sigtab_Hup[row.names(sigtab_Hup) %in% hupAF,]
dim(sigtab_Hup_hupAF)
sigtab_Lup_hupAF <- sigtab_Lup[row.names(sigtab_Lup) %in% hupAF,]
dim(sigtab_Lup_hupAF)

res[which(res$padj < alpha & res$log2FoldChange<0), ]

## upregulated in H
x = tapply(sigtab_Hup$log2FoldChange, sigtab_Hup$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_Hup$Phylum = factor(as.character(sigtab_Hup$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab_Hup$log2FoldChange, sigtab_Hup$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab_Hup$Genus = factor(as.character(sigtab_Hup$Genus), levels=names(x))

# theme_set(theme_bw()) # delete?
# scale_fill_discrete <- function(palname = "Set1", ...) {
#   scale_fill_brewer(palette = palname, ...)
# }

png("logfold_stleaf_Hup_genus031923.png", width=5000, height=3000,res=550)
ggplot(sigtab_Hup, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=5) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  scale_color_manual(values = phylum_col)
dev.off()

#up regulated in L
# Phylum order
x = tapply(sigtab_Lup$log2FoldChange, sigtab_Lup$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_Lup$Phylum = factor(as.character(sigtab_Lup$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab_Lup$log2FoldChange, sigtab_Lup$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab_Lup$Genus = factor(as.character(sigtab_Lup$Genus), levels=names(x))

# 
# theme_set(theme_bw()) #delete?
# scale_fill_discrete <- function(palname = "Set1", ...) {
#   scale_fill_brewer(palette = palname, ...)
# }


png("logfold_stleaf_Lup_genus_031923.png", width=5000, height=3000,res=550)
ggplot(sigtab_Lup, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=5) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  scale_color_manual(values = phylum_col)
dev.off()

saveRDS(res, "res_standard_LleafvsH_oldleaf.RDS")
saveRDS(sigtab_Hup, "stleaf_sigtab_Hup.RDS")
saveRDS(sigtab_Lup, "stleaf_sigtab_Lup.RDS")
saveRDS(sigtab_Hup_hupAF, "stleaf_sigtab_Hup_hupAF.RDS")
saveRDS(sigtab_Lup_hupAF, "stleaf_sigtab_Lup_hupAF.RDS")




#### run correlation test, use abundance data (normalized by proportion) ###
### https://rstudio-pubs-static.s3.amazonaws.com/268156_d3ea37937f4f4469839ab6fa2c483842.html#otus_that_differ_by

# not sure why run this #
colindex <- colnames(otu_table(psf_p_h)) %in% row.names(sigtab_YLup)
otu_psf_p_h_de <- otu_table(psf_p_h)[, colindex]
otu_psf_p_h_upYleaf <- merge_phyloseq(otu_psf_p_h_de, tax_table(psf_p_h), sample_data(psf_p_h))
cor.kendall_bothleaf = cor(otu_table(otu_psf_p_h_upYleaf), sample_data(otu_psf_p_h_upYleaf)$hup_A, method = "kendall")
psf_p_h



####
# Other methods (pearson & spearman will show error message of standard deviation being zero)
cor.kendall = cor(otu_table(psf_p_h), sample_data(psf_p_h)$hup_A, method = "kendall")
which(cor.kendall > 0.4)
ASV_hupA_0.4 <- cor.kendall[which(cor.kendall > 0.4),]
ASV_hupA_0.4_tax <- tax_table(psf_p_h)[names(ASV_hupA_0.4),]

# R_cor.kendall = cor(otu_table(psf_R_h), sample_data(psf_R_h)$hup_A, method = "kendall")
# which(R_cor.kendall > 0.4)
# R_cor.kendall[which(R_cor.kendall > 0.4),]

saveRDS(ASV_hupA_0.4, "ASV_hupA_0.4.RDS")
write.table(ASV_hupA_0.4_tax, "ASV_hupA_0.4_tax.txt", sep = "\t", quote = F, row.names =F)

png("cor_gt_0.4_value.png", width=4000, height=3000,res=400)
par(mfrow = c(2,2))
plot(otu_table(psf_p_h)[,14] ~ sample_data(psf_p_h)$hup_A, main= "Trichothecium roseum τ =0.44", xlab = "relative abundance", ylab = "HupA (μg/g)")
plot(otu_table(psf_p_h)[,10] ~ sample_data(psf_p_h)$hup_A, main="Penicillium citrium τ =0.50", xlab = "relative abundance", ylab = "HupA (μg/g)")
plot(otu_table(psf_p_h)[,1] ~ sample_data(psf_p_h)$hup_A, main="Aspergillus sp. τ = 0.55", xlab = "relative abundance", ylab = "HupA (μg/g)")
plot(otu_table(psf_p_h)[,336] ~ sample_data(psf_p_h)$hup_A, main="Hypocreales sp. τ = 0.44", xlab = "relative abundance", ylab = "HupA (μg/g)")
dev.off()

cor4 <- as.data.frame (otu_table(psf_p_h)[,colnames(otu_table(psf_p_h)) %in% names(ASV_hupA_0.4)])
cor4 <- merge(cor4, data.frame(sample_data(psf_p_h)), by = "row.names")

cor4 <- cor4 %>% 
  pivot_longer(
    cols = starts_with("ASV"),
                       names_to = "ASV", 
                       values_to = "relative_abundance")

png("cor.png", width = 4000, height =3000, res= 500) #Figure terrible
ggplot(data = cor4) +
  geom_point(mapping = aes(x = relative_abundance, y = hup_A, color = ASV)) +
  geom_smooth(
    mapping = aes(x = relative_abundance, y = hup_A, color = ASV),
    show.legend = FALSE
  )
dev.off()


png("cor_gt_0.4_value.png", width=4000, height=3000,res= 550)
par(mfrow = c(2,2))
plot(sample_data(psf_p_h)$hup_A ~ otu_table(psf_p_h)[,colnames(otu_table(psf_p_h)) == names(ASV_hupA_0.4)[1]], main="Trichothecium roseum τ =0.44", xlab = "Relative abundance", ylab = "HupA (μg/g)")
plot(sample_data(psf_p_h)$hup_A ~ otu_table(psf_p_h)[,colnames(otu_table(psf_p_h)) == names(ASV_hupA_0.4)[2]], main="Penicillium citrium τ =0.50", xlab = "Relative abundance", ylab = "HupA (μg/g)")
plot(sample_data(psf_p_h)$hup_A ~ otu_table(psf_p_h)[,colnames(otu_table(psf_p_h)) == names(ASV_hupA_0.4)[3]], main="Aspergillus sp. τ = 0.55", xlab = "Relative abundance", ylab = "HupA (μg/g)")
plot(sample_data(psf_p_h)$hup_A ~ otu_table(psf_p_h)[,colnames(otu_table(psf_p_h)) == names(ASV_hupA_0.4)[4]], main="Hypocreales sp. τ = 0.44", xlab = "Relative abundance", ylab = "HupA (μg/g)")
#plot(otu_table(psf_p_h)[,colnames(otu_table(psf_p_h)) == names(ASV_hupA_0.4)[4]] ~ sample_data(psf_p_h)$hup_A, main="OTU14, Aspergillus sydowii, r= 0.55", xlab = "relative abundance", ylab = "HupA")
dev.off()


# png("cor>0.4_tissue.png", width=4000, height=3000,res=400)
# par(mfrow = c(2,2))
# boxplot(otu_table(psf_p_h)[,14] ~ sample_data(psf_p_h)$tissue, ylab="% Relative abundance", main="OTU14, Aspergillus sydowii, r= 0.55 ")
# boxplot(otu_table(psf_p_h)[,10] ~ sample_data(psf_p_h)$tissue, ylab="% Relative abundance", main="OTU10, Penicillium citrium, r =0.5")
# boxplot(otu_table(psf_p_h)[,1] ~ sample_data(psf_p_h)$tissue, ylab="% Relative abundance", main="OTU1, Trichothecium roseum, r = 0.44")
# boxplot(otu_table(psf_p_h)[,336] ~ sample_data(psf_p_h)$tissue, ylab="% Relative abundance", main="OTU336 Acremonium sp, r = 0.44")
# dev.off()



##### PERMANOVA & alpha diversity test ##
set.seed(1)
# Calculate bray curtis distance matrix
psf_p_h_bray <- phyloseq::distance(psf_p_h, method = "bray")
sampledf <- data.frame(sample_data(psf_p_h))

# Adonis test
# Huperzia only # all tissues
adonis2(psf_p_h_bray ~ tissue, data = sampledf)
pairwise.adonis2(psf_p_h_bray ~ tissue, data = sampledf)

# leafstem of the two species
psf_p_leafstem1 <- psf_p_leafstem
sample_data(psf_p_leafstem1)$tissue[sample_data(psf_p_leafstem1)$tissue %in% c("young_leaf", "standard_leaf")] <- "leaf"
psf_R_leafstem1 <- psf_R_leafstem
sample_data(psf_R_leafstem1)$tissue[sample_data(psf_R_leafstem1)$tissue %in% c("young_leaf", "standard_leaf")] <- "leaf"


psf_p_leafstem1_bray <- phyloseq::distance(psf_p_leafstem1, method = "bray")#check the interaction between species and tissues
adonis2(psf_p_leafstem1_bray ~ tissue*species, data = data.frame(sample_data(psf_p_leafstem1)))

psf_p_leafstem1_l <- subset_samples(psf_p_leafstem1, species == "Lycopodium_complanatum")
psf_p_leafstem1_h <- subset_samples(psf_p_leafstem1, species == "Huperzia_selago")
psf_p_leafstem1_leaf <- subset_samples(psf_p_leafstem1, tissue == "leaf")
psf_p_leafstem1_stem <- subset_samples(psf_p_leafstem1, tissue == "stem")


psf_p_leafstem1_l_bray <- phyloseq::distance(psf_p_leafstem1_l, method = "bray")
psf_p_leafstem1_h_bray <- phyloseq::distance(psf_p_leafstem1_h, method = "bray")
psf_p_leafstem1_leaf_bray <- phyloseq::distance(psf_p_leafstem1_leaf, method = "bray")
psf_p_leafstem1_stem_bray <- phyloseq::distance(psf_p_leafstem1_stem, method = "bray")

adonis2(psf_p_leafstem1_l_bray ~ tissue, data = data.frame(sample_data(psf_p_leafstem1_l)))
adonis2(psf_p_leafstem1_h_bray ~ tissue, data = data.frame(sample_data(psf_p_leafstem1_h)))
adonis2(psf_p_leafstem1_leaf_bray ~ species, data = data.frame(sample_data(psf_p_leafstem1_leaf)))
adonis2(psf_p_leafstem1_stem_bray ~ species, data = data.frame(sample_data(psf_p_leafstem1_stem)))

rich_psf_R_leafstem1 <- data.frame(
  "Observed" = phyloseq::estimate_richness(psf_R_leafstem1, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(psf_R_leafstem1, measures = "Shannon"),
  "tissue" = phyloseq::sample_data(psf_R_leafstem1)$tissue,
  "species" = phyloseq::sample_data(psf_R_leafstem1)$species)

moments::skewness(rich_psf_R_leafstem1$Observed, na.rm = TRUE) #https://bookdown.org/wadetroberts/r-you-ready-for-r/anova.html#packages-needed-for-anova
moments::kurtosis(rich_psf_R_leafstem1$Observed, na.rm = TRUE)
moments::skewness(rich_psf_R_leafstem1$Shannon, na.rm = TRUE) #https://bookdown.org/wadetroberts/r-you-ready-for-r/anova.html#packages-needed-for-anova
moments::kurtosis(rich_psf_R_leafstem1$Shannon, na.rm = TRUE)
shapiro.test(rich_psf_R_leafstem1$Observed)
shapiro.test(rich_psf_R_leafstem1$Shannon)
kruskal.test(Observed ~ tissue, data = rich_psf_R_leafstem1)
kruskal.test(Shannon ~ tissue, data = rich_psf_R_leafstem1)
kruskal.test(Observed ~ species, data = rich_psf_R_leafstem1)
kruskal.test(Shannon ~ species, data = rich_psf_R_leafstem1[rich_psf_R_leafstem1$tissue == "leaf",])
kruskal.test(Observed ~ species, data = rich_psf_R_leafstem1[rich_psf_R_leafstem1$tissue == "leaf",])
kruskal.test(Shannon ~ species, data = rich_psf_R_leafstem1[rich_psf_R_leafstem1$tissue == "stem",])
kruskal.test(Observed ~ species, data = rich_psf_R_leafstem1[rich_psf_R_leafstem1$tissue == "stem",])


rich_psf_R_leafstem1$sample <- row.names(rich_psf_R_leafstem1)
long_rich_psf_R_leafstem1 <- as_tibble(rich_psf_R_leafstem1) %>% 
  pivot_longer(c(`Observed`, `Shannon`), names_to = "Index", values_to = "Values") %>%
  as.data.frame
row.names(long_rich_psf_R_leafstem1) <- long_rich_psf_R_leafstem1$sample


pdf("long_rich_psf_R_leafstem1.pdf", width = 6.3, height = 6)
#png("long_rich_psf_R_leafstem1.png", width = 7300, height = 6000, res = 1200)
p <- ggplot(long_rich_psf_R_leafstem1, aes(x=species, y=Values)) + 
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  facet_wrap(Index ~ tissue, scales = "free_y")+
  theme(axis.text.x = element_text(angle = -60, hjust = 0, vjust=0.1))+
  labs(x="Species", y = "Alpha diversity")
#  geom_jitter(shape=16, position=position_jitter(0.2))
p
dev.off()


##### try lefse for Huperzia different tissues ###
#psf_p_h_taxonly <- psf_p_h
taxonly <- as.data.frame(tibble(as.data.frame(tax_table(psf_p_h))) %>% select(Phylum:Species))
row.names(taxonly) <- row.names(as.data.frame(tax_table(psf_p_h)))
psf_p_h_forlefse <- phyloseq(sample_data(psf_p_h), otu_table(psf_p_h), tax_table(as.matrix(taxonly)))
lef_out<-run_lefse(psf_p_h_forlefse, group = "tissue", norm = "CPM",
                   
                   kw_cutoff = 0.01, lda_cutoff = 2)

lef_out_genus <-run_lefse(psf_p_h_forlefse, group = "tissue", norm = "CPM", taxa_rank = "Genus",
                   
                   kw_cutoff = 0.01, lda_cutoff = 2)


png("lefse_bar_lda4.png", width=4000, height = 14000, res = 550)
plot_ef_bar(lef_out)+
scale_fill_manual(values=tissuecol_all, breaks = c("bulbil", "shoot", "young_leaf", "standard_leaf", "stem","sporangium"))
#color setting doesn't work
dev.off()


png("lefse_bar_genus.png", width=4000, height = 10000, res = 600)
plot_ef_bar(lef_out_genus)+
  scale_fill_manual(values=tissuecol_all, breaks = c("bulbil", "shoot", "young_leaf", "standard_leaf", "stem","sporangium"))
#color setting doesn't work
dev.off()

png("lefse_cladogram.png", width = 50000, height =10000, res = 750)
p <- plot_cladogram(lef_out, color = c(bulbil = "#C77CFF", shoot = "#87CEEB", young_leaf = "#66CC99", standard_leaf = "seagreen4", 
                                   stem = "#FF7F50", sporangium ="#DAA520")) +
  theme(plot.margin = margin(0, 0, 0, 0))
p
dev.off()


##### try lefse for Huperzia vs. Lycopodium ### # later 



######
##############
## Volcano plot ## doesn't seem working ## https://lashlock.github.io/compbio/R_presentation.html
####################33
png("test", width=5000, height=3000,res=550)
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()
##


####################### SAVE dataset ##
#######################

saveRDS(ps, "ps.rds")
saveRDS(psf, "psf.rds")
saveRDS(lu_psf, "lu_psf.rds")
saveRDS(psf_p, "psf_p.rds") #dataset removed singleton, normalized by proportion
saveRDS(psf_R, "psf_R.rds") #dataset kept singleton, normalized by rarefaction


## detect differentially present OTUs
psf_p <- readRDS("psf_p.rds")
psf_p_leaf <- subset_samples(psf_p, tissue == "leaf")

test <- simper(otu_table(psf_p_leaf), sample_data(psf_p_leaf)$species, permutations=100)
str(test)


