#/usr/bin/Rscript

start_time <- Sys.time()

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(tidytree))

suppressPackageStartupMessages(library(treeio, include.only="isTip"))
suppressPackageStartupMessages(library(phytools, include.only="midpoint.root"))
suppressPackageStartupMessages(library(ape, include.only="write.nexus"))
suppressPackageStartupMessages(library(gplots, include.only="col2hex"))


###### PARAMETERS AND SOURCE FUNCTIONS

root = 0
dev=0
#dev=1
if(dev) {
	source("bin/functions.R")
	arguments <- list()
}

run <- commandArgs(trailingOnly=FALSE)
file <- grep("--file=",run)
filename <- sub("--file=","",run[file])
path <- dirname(filename)
functionsfile <- paste0(path,"/functions.R")
source(functionsfile)

end_time1 <- Sys.time()
secs <- as.numeric(end_time1-start_time) %>% round(digits=2)
catyellow(paste0("Loading time: ",secs, " seconds"))


##### DATA FILES

arguments <- readArgs(run)

treefile <- arguments$tree
taxafile <- arguments$groups
colorfile <- arguments$color
groupFromNames <- arguments$groupFromName
separator <- arguments$sep
suffix <- arguments$suffix
ignorePrefix <- arguments$ignorePrefix
threads <- arguments$threads
outfileCol <- arguments$outfileCol
outfileColNxs <- arguments$outfileColNxs
outfileUncol <- arguments$outfileUncol
outfileUncolNxs <- arguments$outfileUncolNxs
midpoint <- arguments$midpoint
widthParam <- ifelse(length(arguments$width) > 0, arguments$width, 7)
widthParamCol <- as.numeric(arguments$widthCol)
widthParamUncol <- as.numeric(arguments$widthUncol)


###### READ DATA

catyellow("Reading tree file...")
tree <- read.tree(treefile)
if (length(midpoint) > 0) {
	tree <- midpoint.root(tree)
}

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

col_groups <- tibble(group=character(),col=character())
if (length(colorfile) > 0) {
	catyellow("Reading taxon colors file...")
	col_groups <- read_tsv(
			colorfile,
			col_select = c(1,2),
			col_names = c("group", "col"),
			show_col_types = FALSE
		   )
	cat("\n")
	col_groups <- col_groups %>%
		add_row(group="notfound", col="gray15")
}

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
outtreeCol <- objects[[3]]
outtreeUncol <- objects[[4]]

write.nexus(outtreeCol, file=outfileColNxs)
write.nexus(outtreeUncol, file=outfileUncolNxs)


###### GROUP NODES REPRESENTING MONOPHYLETIC GROUPS (REQUIRED FOR COLLAPSING)

groupedTree <- groupClade(tree, monoNodes$node)
nNodes <- groupedTree$Nnode

###### CREATE BASE TREE, RESCALE MONOPHYLETIC GROUPS AND PRINT THEM

end_time2 <- Sys.time()
secs <- as.numeric(end_time2-end_time1) %>% round(digits=2)
catyellow(paste0("Processing time: ",secs, " seconds"))

cat("\n")
catyellow("Preparing collapsed tree...")
cat("\n")

p <- ggtree(groupedTree) + xlim(NA,x_limit)
p <- collapse_treeplot(p, monoNodes, nNodes)
if(length(color_vector_groups)>0) {
	p <- color_branches(p, groupedOTUs, color_vector_groups)
}
p <- draw_support_values(p)
p <- print_tiplabels(p)
p <- draw_root(tree, p, root)
p <- p + theme_tree2(legend.position="bottom") + guides(color="none")
#p <- p + theme_tree2() + guides(color="none")


###### PRINT
catyellow("Plotting collapsed tree...")

heightCollapsed <- calculateHeightCollapsed(tree, monoNodes)
widthCol <- ifelse(length(widthParamCol) > 0, widthParamCol, widthParam)


cairo_pdf(outfileCol, family="Helvetica",height=heightCollapsed,width=widthCol)
p
invisible(dev.off())


#####################################
#####################################
##   Tree with uncollapsed groups  ##
#####################################
#####################################


catyellow("Plotting uncollapsed tree...")

q <- ggtree(tree) + xlim(NA,x_limit)
if(length(color_vector_groups)>0) {
	q <- color_branches(q, groupedOTUs, color_vector_groups)
}
q <- print_support_values(q)
q <- draw_root(tree, q, root)
#q <- q + theme_tree2()
q <- q + theme_tree2(legend.position="bottom")

heightUncollapsed <- calculateHeightUncollapsed(tree)
widthUncol <- ifelse(length(widthParamUncol) > 0, widthParamUncol, widthParam)

pdf(outfileUncol,height=heightUncollapsed,width=widthUncol)
q
invisible(dev.off())

end_time3 <- Sys.time()
secs <- as.numeric(end_time3-end_time2) %>% round(digits=2)
catyellow(paste0("Plotting time: ",secs, " seconds"))


end_time <- Sys.time()
secs <- as.numeric(end_time-start_time,units="secs") %>% round(digits=2)
mins <- as.numeric(end_time-start_time,units="mins") %>% round(digits=2)

cat("\n")
catyellow(paste0("Total time: ",secs, " seconds (", mins, " mins)"))
catyellow("Enjoy your tree! =)")
cat("\n")
