# Tutorial Setup File for all Dependencies
# =Modified from=
# https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/

#NOTE: When running code, besure to open an interactive R terminal session and run these block-by-block for best results. 
# Highly recommend highlighting whole block and running by CTRL+ENTER (relevant in VSCode)

#%# Slide 3 #%#
## Installing packages
.cran_packages <- c("tidyverse", "cowplot", "picante", "vegan", "HMP", "dendextend", "rms", "devtools")
.bioc_packages <- c("phyloseq", "DESeq2", "microbiome", "metagenomeSeq", "ALDEx2")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(.bioc_packages, version = "3.12") #Different from the website's "3.9" since I'm using R 4.0.5
devtools::install_github("adw96/breakaway")
devtools::install_github(repo = "UVic-omics/selbal")

## Loading Libraries
library(tidyverse); packageVersion("tidyverse")    
library(phyloseq); packageVersion("phyloseq")   
library(DESeq2); packageVersion("DESeq2")
library(microbiome); packageVersion("microbiome") 
library(vegan); packageVersion("vegan")  
library(picante); packageVersion("picante") 
library(ALDEx2); packageVersion("ALDEx2") 
library(metagenomeSeq); packageVersion("metagenomeSeq") 
library(HMP); packageVersion("HMP")  
library(dendextend); packageVersion("dendextend")  
library(selbal); packageVersion("selbal")  
library(rms); packageVersion("rms")
library(breakaway); packageVersion("breakaway") 


#%# Slide 6 #%#
#Get data
(ps <- readRDS("ps_giloteaux_2016.rds"))
sort(phyloseq::sample_sums(ps)) 


#%# Slide 7 #%#
(ps <- phyloseq::subset_samples(ps, phyloseq::sample_sums(ps) > 5000)) 

(ps <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps) > 0, ps)) 


#%# Slide 8 #%#
phyloseq::sample_data(ps)$Status <- ifelse(phyloseq::sample_data(ps)$Subject == "Patient", "Chronic Fatigue", "Control")
phyloseq::sample_data(ps)$Status <- factor(phyloseq::sample_data(ps)$Status, levels = c("Control", "Chronic Fatigue"))
ps %>% 
    sample_data %>% 
    dplyr::count(Status)

#%# Slide 9 #%#
#Visualization of the data
table(phyloseq::tax_table(ps)[, "Phylum"])
ps_rel_abund = phyloseq::transform_sample_counts(ps, function(x){x / sum(x)})

# Operational Taxonomic Unit (OTU)
phyloseq::otu_table(ps)[1:5, 1:5]
phyloseq::otu_table(ps_rel_abund)[1:5, 1:5]


#%# Slide 10 #%#
phyloseq::plot_bar(ps_rel_abund, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Status, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


#%# Slide 11 #%#
ps_phylum <- phyloseq::tax_glom(ps, "Phylum")
phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
phyloseq::otu_table(ps_phylum)[1:5, 1:5]


#%# Slide 12 #%#
#Melt then plot
phyloseq::psmelt(ps_phylum) %>%
ggplot(data = ., aes(x = Status, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")


#%# Slide 13 #%#
#Subset groups
controls <- phyloseq::subset_samples(ps_phylum, Status == "Control")
cf <- phyloseq::subset_samples(ps_phylum, Status == "Chronic Fatigue")
#Output OTU tables
control_otu <- data.frame(phyloseq::otu_table(controls))
cf_otu <- data.frame(phyloseq::otu_table(cf))
#Group rare phyla
control_otu <- control_otu %>%
  t(.) %>%
  as.data.frame(.) %>%
  mutate(Other = Cyanobacteria + Euryarchaeota + Tenericutes + Verrucomicrobia + Fusobacteria) %>%
  dplyr::select(-Cyanobacteria, -Euryarchaeota, -Tenericutes, -Verrucomicrobia, -Fusobacteria)
cf_otu <- cf_otu %>%
  t(.) %>%
  as.data.frame(.) %>%
  mutate(Other = Cyanobacteria + Euryarchaeota + Tenericutes + Verrucomicrobia + Fusobacteria) %>%
  dplyr::select(-Cyanobacteria, -Euryarchaeota, -Tenericutes, -Verrucomicrobia, -Fusobacteria)
#HMP test
group_data <- list(control_otu, cf_otu)
(xdc <- HMP::Xdc.sevsample(group_data)) 


#%# Slide 14 #%#
#Hierarchical Clustering
ps_rel_otu <- data.frame(phyloseq::otu_table(ps_rel_abund))
ps_rel_otu <- t(ps_rel_otu)
bc_dist <- vegan::vegdist(ps_rel_otu, method = "bray")
as.matrix(bc_dist)[1:5, 1:5]


#%# Slide 15 #%#
#Save as dendrogram
ward <- as.dendrogram(hclust(bc_dist, method = "ward.D2"))
#Provide color codes
meta <- data.frame(phyloseq::sample_data(ps_rel_abund))
colorCode <- c(Control = "red", `Chronic Fatigue` = "blue")
labels_colors(ward) <- colorCode[meta$Status][order.dendrogram(ward)]
#Plot
plot(ward)

#%# Slide 16 #%#
#Alpha-Diversity
ggplot(data = data.frame("total_reads" =  phyloseq::sample_sums(ps),
                         "observed" = phyloseq::estimate_richness(ps, measures = "Observed")[, 1]),
       aes(x = total_reads, y = observed)) +
  geom_point() +
  geom_smooth(method="lm", se = FALSE) +
  labs(x = "\nTotal Reads", y = "Observed Richness\n")


#%# Slide 17 #%#
#Subsample reads
(ps_rare <- phyloseq::rarefy_even_depth(ps, rngseed = 123, replace = FALSE)) 
head(phyloseq::sample_sums(ps_rare))
adiv <- data.frame(
  "Observed" = phyloseq::estimate_richness(ps_rare, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps_rare, measures = "Shannon"),
  #Add include.root=FALSE
  "PD" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(ps_rare)))), tree = phyloseq::phy_tree(ps_rare), include.root=FALSE)[, 1],
  "Status" = phyloseq::sample_data(ps_rare)$Status)
head(adiv)

#%# Slide 18 #%#
#Plot adiv measures
adiv %>%
  gather(key = metric, value = value, c("Observed", "Shannon", "PD")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "PD"))) %>%
  ggplot(aes(x = Status, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Status), height = 0, width = .2) +
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  theme(legend.position="none")

#Summarize
adiv %>%
  group_by(Status) %>%
  dplyr::summarise(median_observed = median(Observed),
            median_shannon = median(Shannon),
            median_pd = median(PD))

#%# Slide 19 #%#
#Wilcoxon test of location
wilcox.test(Observed ~ Status, data = adiv, exact = FALSE, conf.int = TRUE)

wilcox.test(Shannon ~ Status, data = adiv, conf.int = TRUE)    

wilcox.test(PD ~ Status, data = adiv, conf.int = TRUE)


#%# Slide 20 #%#
#Obtain breakaway estimates
ba_adiv <- breakaway::breakaway(ps)
ba_adiv[1]
#Examine models
summary(ba_adiv) %>%
  add_column("SampleNames" = ps %>% otu_table %>% sample_names)  
#Test for group differnce
bt <- breakaway::betta(summary(ba_adiv)$estimate,
                       summary(ba_adiv)$error,
                       make_design_matrix(ps, "Status"))
bt$table 


#%# Slide 21 #%#
#Plot estimates
plot(ba_adiv, ps, color = "Status")     


#%# Slide 22
#Beta-Diversity
#CLR transform
(ps_clr <- microbiome::transform(ps, "clr"))  

phyloseq::otu_table(ps)[1:5, 1:5]

phyloseq::otu_table(ps_clr)[1:5, 1:5]


#%# Slide 23 #%#
#PCA via phyloseq
ord_clr <- phyloseq::ordinate(ps_clr, "RDA")
#Plot scree plot
phyloseq::plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")


#%# Slide 24 #%#
#Examine eigenvalues and % prop. variance explained
head(ord_clr$CA$eig)     

sapply(ord_clr$CA$eig[1:5], function(x) x / sum(ord_clr$CA$eig))    


#%# Slide 25 #%#
#Scale axes and plot ordination
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
phyloseq::plot_ordination(ps, ord_clr, type="samples", color="Status") + 
  geom_point(size = 2) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = Status), linetype = 2)


#%# Slide 26 #%#
#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(ps_clr, method = "euclidean") 
#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$Status)


#%# Slide 27 #%#
#Dispersion test and plot
dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(ps_clr)$Status)
dispr

plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")


#%# Slide 28 #%#
boxplot(dispr, main = "", xlab = "")

permutest(dispr)


#%# Slide 29 #%#
#Generate distances
ord_unifrac <- ordinate(ps_rare, method = "PCoA", distance = "wunifrac") 
ord_unifrac_un <- ordinate(ps_rare, method = "PCoA", distance = "unifrac")   
#Plot ordinations
a <- plot_ordination(ps_rare, ord_unifrac, color = "Status") + geom_point(size = 2)
b <- plot_ordination(ps_rare, ord_unifrac_un, color = "Status") + geom_point(size = 2)
cowplot::plot_grid(a, b, nrow = 1, ncol = 2, scale = .9, labels = c("Weighted", "Unweighted"))


#%# Slide 30 #%#
#Differential Abundance Testing
#Generate data.frame with OTUs and metadata
ps_wilcox <- data.frame(t(data.frame(phyloseq::otu_table(ps_clr))))
ps_wilcox$Status <- phyloseq::sample_data(ps_clr)$Status
#Define functions to pass to map
wilcox_model <- function(df){
  wilcox.test(abund ~ Status, data = df)
}
wilcox_pval <- function(df){
  wilcox.test(abund ~ Status, data = df)$p.value
}
#Create nested data frames by OTU and loop over each using map 
wilcox_results <- ps_wilcox %>%
  gather(key = OTU, value = abund, -Status) %>%
  group_by(OTU) %>%
  nest() %>%
  mutate(wilcox_test = map(data, wilcox_model),
         p_value = map(data, wilcox_pval))                       

#%# Slide 31 #%#
#Show results
head(wilcox_results)

head(wilcox_results$data[[1]])

wilcox_results$wilcox_test[[1]]

wilcox_results$p_value[[1]]


#%# Slide 32 #%#
#Unnesting
wilcox_results <- wilcox_results %>%
  dplyr::select(OTU, p_value) %>%
  unnest()

head(wilcox_results)

#Adding taxonomic labels
taxa_info <- data.frame(tax_table(ps_clr))
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")
#Computing FDR corrected p-values
wilcox_results <- wilcox_results %>%
  full_join(taxa_info) %>%
  arrange(p_value) %>%
  mutate(BH_FDR = p.adjust(p_value, "BH")) %>%
  filter(BH_FDR < 0.05) %>%
  dplyr::select(OTU, p_value, BH_FDR, everything())

#Printing results
print.data.frame(wilcox_results)  


#%# Slide 33 #%#
#Run ALDEx2 (runs for a little bit)
aldex2_da <- ALDEx2::aldex(data.frame(phyloseq::otu_table(ps)), phyloseq::sample_data(ps)$Status, test="t", effect = TRUE, denom="iqlr")

#Plot effect sizes
ALDEx2::aldex.plot(aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)


#%# Slide 34 #%#
#Clean up presentation
sig_aldex2 <- aldex2_da %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, wi.ep, wi.eBH)
sig_aldex2 <- left_join(sig_aldex2, taxa_info)

sig_aldex2


#%# Slide 35 #%#
#Prediction
#Generate data.frame
clr_pcs <- data.frame(
  "pc1" = ord_clr$CA$u[,1],
  "pc2" = ord_clr$CA$u[,2],
  "pc3" = ord_clr$CA$u[,3],
  "Status" = phyloseq::sample_data(ps_clr)$Status
)
clr_pcs$Status_num <- ifelse(clr_pcs$Status == "Control", 0, 1)
head(clr_pcs)

#Specify a datadist object (for rms)
dd <- datadist(clr_pcs)
options(datadist = "dd")

#Plot the unconditional associations
a <- ggplot(clr_pcs, aes(x = pc1, y = Status_num)) +
  Hmisc::histSpikeg(Status_num ~ pc1, lowess = TRUE, data = clr_pcs) +
  labs(x = "\nPC1", y = "Pr(Chronic Fatigue)\n")
b <- ggplot(clr_pcs, aes(x = pc2, y = Status_num)) +
  Hmisc::histSpikeg(Status_num ~ pc2, lowess = TRUE, data = clr_pcs) +
  labs(x = "\nPC2", y = "Pr(Chronic Fatigue)\n")
c <- ggplot(clr_pcs, aes(x = pc3, y = Status_num)) +
  Hmisc::histSpikeg(Status_num ~ pc3, lowess = TRUE, data = clr_pcs) +
  labs(x = "\nPC3", y = "Pr(Chronic Fatigue)\n")
cowplot::plot_grid(a, b, c, nrow = 2, ncol = 2, scale = .9, labels = "AUTO")


#%# Slide 36 #%#
#Fit full model with splines (3 knots each)
m1 <- rms::lrm(Status_num ~ rcs(pc1, 3) + rcs(pc2, 3) + rcs(pc3, 3), data = clr_pcs, x = TRUE, y = TRUE)

#Grid search for penalties
pentrace(m1, list(simple = c(0, 1, 2), nonlinear = c(0, 100, 200)))

pen_m1 <- update(m1, penalty = list(simple = 1, nonlinear = 200))
pen_m1


#%# Slide 37 #%#
#Plot log odds
ggplot(Predict(pen_m1))


#%# Slide 38 #%#
#Obtain optimism corrected estimates
(val <- rms::validate(pen_m1))

#Compute corrected c-statistic
(c_opt_corr <- 0.5 * (val[1, 5] + 1))


#%# Slide 39 #%#
#Plot calibration
cal <- rms::calibrate(pen_m1, B = 200)
plot(cal)

#Output pred. probs
head(predict(pen_m1, type ="fitted"))

#%# Slide 40 #%#
#Selbel alternative
#Agglomerate taxa
(ps_family <- phyloseq::tax_glom(ps, "Family"))

#Windows Defender Might Complain a bit
phyloseq::taxa_names(ps_family) <- phyloseq::tax_table(ps_family)[, "Family"]
#Run selbal
cv_sebal <- selbal::selbal.cv(x = data.frame(t(data.frame(phyloseq::otu_table(ps_family)))), 
                              y = phyloseq::sample_data(ps_family)$Status, 
                             n.fold = 5, n.iter = 1) 


#%# Slide 41 #%#
#plot/print results
cv_sebal$accuracy.nvar


#%#% Slide 42 #%#
plot.new()
grid.draw(cv_sebal$global.plot)