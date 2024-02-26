# Setup ####
library(tidyverse)
library(phyloseq)
library(rbiom)
library(ape)
library(TreeTools)
library(vegan)
set.seed(666)

# Load data ####

# phylogeny (newick format)
tree <- ape::read.tree("./data/tree.nwk")
# asv table (from DADA2 | sample IDs stored as first column named "sample_id")
asv <- read.csv("./data/asv_table.csv",row.names = "sample_id")
# metadata
meta <- read_csv("./data/metadata.csv")

# re-root tree (picking a random outgroup for demo)
tree <- TreeTools::RootTree(tree, outgroupTips = "355132")


# Check data ####
# check that tip labels match order of taxa in asv table
# if not, rearrange asv rows to match
if(!identical(row.names(asv),tree$tip.label)){
  asv <- asv[tree$tip.label,]
}


# UniFrac distance ####

# normalize counts in each sample
asv_relabund <- apply(as.matrix(asv),2,function(x){x/sum(x)})

# rbiom::unifrac() assumes that taxa are rows and samples are columns
dist_wu <- rbiom::unifrac(biom=asv_relabund,weighted = TRUE,tree=tree)

# Ordinate ####
# NMDS with unifrac distance
nmds <- vegan::monoMDS(dist_wu,model = "global")

# add ordination points to metadata for modeling/plotting
meta$mds1 <- nmds$points[,1]
meta$mds2 <- nmds$points[,2]


# plot ####
meta %>% 
  ggplot(aes(x=mds1,y=mds2,color=SampleType)) +
  geom_point()

# Try in phyloseq ####
met <- sample_data(meta)
sample_names(met) <- met$X.SampleID

ps <- phyloseq(met,
               otu_table(asv,taxa_are_rows = TRUE),
               phy_tree(tree))

ord <- ps %>% ordinate(method="NMDS",distance = "unifrac")
plot_ordination(ps, ord, color="SampleType")

# honestly, phyloseq version looks nicer...