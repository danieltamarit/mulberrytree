#/usr/bin/Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(ggtree))
suppressMessages(library(tidytree))


###### PARAMETERS

root = 0


##### DATA FILES

run <- commandArgs(trailingOnly=FALSE)
file <- grep("--file=",run)
filename <- sub("--file=","",run[file])
path <- dirname(filename)
functionsfile <- paste0(path,"/bin/functions.R")
source(functionsfile)


args <- readArgs(run)

treefile <- args[[1]]
taxafile <- args[[2]]
colorfile <- args[[3]]
groupFromNames <- args[[4]]

###### READ DATA

catyellow("Reading tree file...")
tree <- read.tree(treefile)

taxaFromFile <- tibble()
if ((length(taxafile) > 0) && (file.exists(taxafile))) {
	catyellow("Reading taxon names file...")
	taxaFromFile <- read_tsv(taxafile,col_names=c("name","group"),show_col_types = FALSE)
}
taxa <- tibble()
if (length(taxaFromFile) > 0) {
	taxa <- taxaFromFile
}
if (groupFromNames == "yes") {
	taxa <- taxaFromNames(tree,taxa)
}


catyellow("Reading taxon colors file...")
col_groups <- read_tsv(
		colorfile,
		col_names = c("group", "col"),
		show_col_types = FALSE
	   )
cat("\n")
catyellow("Analysing tree...")

###### PROCESS TAXON NAMES


#full_names <- tibble(leaves=tree$tip.label, name=sub("_[0-9]+G$", "", tree$tip.label, perl=TRUE))
full_names <- tibble(leaves=tree$tip.label, name=tree$tip.label)
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
