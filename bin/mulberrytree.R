#/usr/bin/Rscript

start_time <- Sys.time()

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
functionsfile <- paste0(path,"/functions.R")
source(functionsfile)

arguments <- readArgs(run)

treefile <- arguments$tree
taxafile <- arguments$groups
colorfile <- arguments$color
groupFromNames <- arguments$groupFromName
separator <- arguments$sep
suffix <- arguments$suffix
threads <- arguments$threads
outfile <- arguments$outfile

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
if ((length(groupFromNames) > 0) && (groupFromNames == "yes")) {
	taxa <- taxaFromNames(tree,taxa,separator,suffix)
}


catyellow("Reading taxon colors file...")
col_groups <- read_tsv(
		colorfile,
		col_names = c("group", "col"),
		show_col_types = FALSE
	   )
cat("\n")


if ((length(threads)==0) || (threads == "")) {
	threads <- 1
} else {
	suppressMessages(library(parallel, include.only="mclapply"))
}

###### PROCESS TAXON NAMES
catyellow("Analysing tree...")

leaves <- getLeafNames(tree, taxa, suffix)

###### PROCESS COLORS

color_vector_groups <- col_groups$col
names(color_vector_groups) <- col_groups$group

###### NORMALISE X AXIS OF THE PLOT
x_limit <- calculate_x_limit(tree)



##### GATHER MONOPHYLY DATA TO COLLAPSE MONOPHYLETIC GROUPS

objects <- monophyletic_subgroups(tree, leaves, col_groups, threads)
monoNodes <- objects[[1]]
groupedOTUs <- objects[[2]]


###### GROUP NODES REPRESENTING MONOPHYLETIC GROUPS (REQUIRED FOR COLLAPSING)

groupedTree <- groupClade(tree, monoNodes$node)
nNodes <- groupedTree$Nnode

###### CREATE BASE TREE, RESCALE MONOPHYLETIC GROUPS AND PRINT THEM

cat("\n")
catyellow("Preparing collapsed tree...")
cat("\n")

p <- ggtree(groupedTree) + xlim(NA,x_limit)
p <- collapse_tree(p, monoNodes, nNodes)
p <- color_branches(p, groupedOTUs, color_vector_groups)
p <- draw_support_values(p)
p <- print_tiplabels(p)
p <- draw_root(tree, p, root)
p <- p + theme_tree2() + guides(color="none")


###### PRINT
catyellow("Plotting collapsed tree...")

heightCollapsed <- calculateHeightCollapsed(tree, monoNodes)

outfile_collapsed = paste0(outfile,"_collapsed.pdf")
cairo_pdf(outfile_collapsed, family="Liberation Sans",height=heightCollapsed)
p
invisible(dev.off())


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

heightUncollapsed <- calculateHeightUncollapsed(tree)

outfile_uncollapsed = paste0(outfile,"_uncollapsed.pdf")
pdf(outfile_uncollapsed,height=heightUncollapsed)
q
invisible(dev.off())

end_time <- Sys.time()
secs <- as.numeric(end_time-start_time) %>% round(digits=2)
mins <- as.numeric(end_time-start_time,units="mins") %>% round(digits=2)

cat("\n")
catyellow(paste0("Total time: ",secs, " seconds (", mins, " mins)"))
catyellow("Enjoy your tree! =)")
cat("\n")
