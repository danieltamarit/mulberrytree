#/usr/bin/Rscript

run <- commandArgs(trailingOnly=FALSE)
file <- grep("--file=",run)
filename <- sub("--file=","",run[file])
path <- dirname(filename)
functionsfile <- paste0(path,"/bin/functions.R")
source(functionsfile)

suppressMessages(library(tidyverse))
suppressMessages(library(ggtree))
suppressMessages(library(tidytree))


###### PARAMETERS

root = 0



##### DATA FILES

#treefile <- "NM-A.fasta.PMSF_LGC60G4F.treefile.renamed.rooted"
#taxafile <- "/local/one/dtamarit/asgard/phylogenomics/trees_laura/All_results_organised_newData/z_classif_all.tsv"
#colorfile <- "/local/one/dtamarit/asgard/phylogenomics/trees_laura/All_results_organised_newData/label_coloring/z_color_groups.tsv"

args <- readArgs(run)

treefile <- args[[1]]
taxafile <- args[[2]]
colorfile <- args[[3]]

#"/local/one/dtamarit/asgard/phylogenomics/alpaca_phylogenomics/220824_testtrees/results/NM_testNew_noDPANN_noKor_221t_20171s.fa.SR6.PMSF_GTRG4C60SR6F_2.treefile",

# args <- c(
# 	"/local/one/dtamarit/asgard/phylogenomics/alpaca_phylogenomics/220824_testtrees/results/NM_testNew_noDPANN_noKor_195t_19277s.fa.PMSF_LGR4C60F.treefile",
# 	"/local/one/dtamarit/asgard/phylogenomics/alpaca_phylogenomics/220824_testtrees/z_classif_all.tsv",
# 	"/local/one/dtamarit/asgard/phylogenomics/alpaca_phylogenomics/220824_testtrees/z_color_groups.tsv"
# )
# treefile <- args[1]
# taxafile <- args[2]
# colorfile <- args[3]

###### READ DATA

catyellow("Reading tree file...")
tree <- read.tree(treefile)

catyellow("Reading taxon names file...")
taxa <- read_tsv(taxafile,col_names=c("name","group"),show_col_types = FALSE)

catyellow("Reading taxon colors file...")
col_groups <- read_tsv(
		colorfile,
		col_names = c("group", "col"),
		show_col_types = FALSE
	   )

cat("\n")
catyellow("Analysing tree...")

###### PROCESS TAXON NAMES


full_names <- tibble(leaves=tree$tip.label, name=sub("_[0-9]+G$", "", tree$tip.label, perl=TRUE))
leaves <- inner_join(full_names, taxa, by='name')

###### PROCESS COLORS

color_vector_groups <- col_groups$col
names(color_vector_groups) <- col_groups$group

###### NORMALISE X AXIS OF THE PLOT
x_limit <- calculate_x_limit(tree)



##### GATHER MONOPHYLY DATA TO COLLAPSE MONOPHYLETIC GROUPS

objects <- monophyletic_subgroups(tree, leaves, col_groups)
monoNodes <- objects[[1]]
groupedOTUs <- objects[[2]]


###### GROUP NODES REPRESENTING MONOPHYLETIC GROUPS (REQUIRED FOR COLLAPSING)

groupedTree <- groupClade(tree, monoNodes$node)
nNodes <- groupedTree$Nnode

###### BASE TREE, RESCALE MONOPHYLETIC GROUPS AND PRINT THEM

cat("\n")
catyellow("Plotting collapsed tree...")

p <- ggtree(groupedTree) + xlim(NA,x_limit)
p <- collapse_tree(p, monoNodes, nNodes)
p <- color_branches(p, groupedOTUs, color_vector_groups)
p <- draw_support_values(p)
p <- print_tiplabels(p)
p <- draw_root(tree, p, root)
p <- p + theme_tree2() + guides(color="none")


###### PRINT

outfile_collapsed = paste0(treefile,"_collapsed.pdf")
cairo_pdf(outfile_collapsed, family="Liberation Sans")
p
invisible(dev.off())






#quit(save="no")



#####################################
#####################################
##   Tree with uncollapsed groups  ##
#####################################
#####################################

catyellow("Plotting uncollapsed tree...")

q <- ggtree(tree) + xlim(NA,x_limit)
q <- color_branches(q, groupedOTUs, color_vector_groups)
q <- print_support_values(q)
q <- draw_root(tree, q, root)
q <- q + theme_tree2()

outfile_uncollapsed = paste0(treefile,"_uncollapsed.pdf")
pdf(outfile_uncollapsed,height=12)
q
invisible(dev.off())



######################################
######################################
##    Tree with collapsed groups    ##
# (only if single monophyletic group)#
######################################
######################################

##### GATHER MONOPHYLY DATA TO COLLAPSE MONOPHYLETIC GROUPS

#objects <- monophyleticGroups(tree, leaves, col_groups)
#monoNodes <- objects[[1]]
#groupedOTUs <- objects[[2]]

###### GROUP NODES REPRESENTING MONOPHYLETIC GROUPS (REQUIRED FOR COLLAPSING)

#groupedTree <- groupClade(tree, monoNodes$node)
#nNodes <- groupedTree$Nnode

###### BASE TREE, RESCALE MONOPHYLETIC GROUPS AND PRINT THEM

#p <- ggtree(groupedTree) + xlim(NA,x_limit)
#p <- rescaleTree(p, monoNodes, nNodes)
#p <- color_branches_in_plot(p, groupedOTUs, color_vector_groups)
#p <- draw_support_values(p)

###### PRINT

#outfile_collapsed = paste0(treefile,"_collapsed.pdf")
#pdf(outfile_collapsed)
#p + theme_tree2()
#invisible(dev.off())



######################################
######################################
##    Tree with collapsed groups    ##
# (only if single monophyletic group)#
######################################
######################################

##### GATHER MONOPHYLY DATA TO COLLAPSE MONOPHYLETIC GROUPS
