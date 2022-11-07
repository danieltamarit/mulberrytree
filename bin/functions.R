#/usr/bin/Rscript

##### ARGUMENTS
readArgs <- function(run) {
#   argtext <- sub(".*--args ",runline)
#   args <- unlist(strsplit(argtext," "))
   args <- run[grep("^--",run, invert=TRUE)]
   args <- args[-1]

   do_not_match <- grep("^tree=|^groups=|^colors=|^groupFromName=|^sep=|^suffix=|^threads=", args, invert=TRUE)
   if (length(do_not_match) > 0) {
      catyellow("#####################################################\n")
      catyellow(
         paste(
             "\nCould not understand the following arguments:\n", paste(args[do_not_match],collapse="\n"),
             "\n\nPlease check your command.\n\nYou can run this script without arguments to read the usage information:\n Rscript mulberrytree.R\n\n", sep=""
          )
       )
       catyellow("#####################################################\n")
       stop()
   }

   tree_arg <- grep("^tree=", args)
   treefile <- unlist(strsplit( args[ tree_arg ], "="))[2]
   if (! file.exists(treefile)) {
      stop("Did not find tree file!")
   }

   group_arg <- grep("^groups=", args)
   groupfile <- unlist(strsplit( args[ group_arg ], "="))[2]
   if ((length(groupfile) > 0) && (! file.exists(groupfile))) {
      stop(paste0("Group file provided but not found: ", groupfile))
   }

   color_arg <- grep("^colors=", args)
   colorfile <- unlist(strsplit( args[ color_arg ], "="))[2]
   if ((length(colorfile) > 0) && (! file.exists(colorfile))) {
      stop(paste0("Color file provided but not found: ", colorfile))
   }

   groupFromName_arg <- grep("^groupFromName=", args)
   groupFromName <- unlist(strsplit( args[ groupFromName_arg ], "="))[2]

   separator_arg <- grep("^sep=", args)
   separator <- unlist(strsplit( args[ separator_arg ], "="))[2]
   if (! length(separator) > 0) {
      separator <- "|"
   }

   suffix_arg <- grep("^suffix=", args)
   suffix <- unlist(strsplit( args[ suffix_arg ], "="))[2]

   threads_arg <- grep("^threads=", args)
   threads <- unlist(strsplit( args[ threads_arg ], "="))[2]

   catyellow("\nmulberrytree will use the following information:")
   catyellow("---------------------------------------")
   catyellow("-Tree file: ")
   catcyan(treefile)
   if ((length(groupfile) > 0) && (file.exists(groupfile))) {
      catyellow("-Group file: ")
      catcyan(groupfile)
   }
   if ((length(colorfile) > 0) && (file.exists(colorfile))) {
      catyellow("-Color file: ")
      catcyan(colorfile)
   }
   if (length(groupFromName) > 0) {
      if (groupFromName == "yes") {
         if (length(groupfile) > 0) {
            if (file.exists(groupfile)) {
               catyellow(paste0(
                           "-Group interpretation from leaf names unless clashing with \"",
                           groupfile,
                           "\".\n-Separator:"
                        ))
               cat(paste0("\"",separator,"\"\n"))
            }
            if (! file.exists(groupfile)) {
               catyellow("-Group interpretation from leaf names.\n-Separator:")
               cat(paste0("\"",separator,"\"\n"))
            }
         }
         if (length(groupfile) == 0) {
            catyellow("-Group interpretation from leaf names.\n-Separator:")
            cat(paste0("\"",separator,"\"\n"))
         }
      }
   }
   if (length(suffix) > 0) {
      catyellow("-Suffix to ignore in tree leaf names:")
      cat(paste0("\"",suffix,"\"\n"))
   }
   if (length(threads) > 0) {
      catyellow("-Number of threads:")
      cat(paste0("\"",threads,"\"\n"))
   }

   catyellow("--------------------------------------\n")



   paramlist <- list(
      tree=treefile,
      groups=groupfile,
      color=colorfile,
      groupFromName=groupFromName,
      sep=separator,
      suffix=suffix,
      threads=threads
   )
   return(paramlist)
}


##### PRINTING FUNCTIONS
catgreen <- function(text) {
   col_start <- "\033[0;0;32m"
   col_end <- "\033[0m"
   cat(paste0(col_start, text, col_end, "\n"))
}
catyellow <- function(text) {
   col_start <- "\033[0;0;33m"
   col_end <- "\033[0m"
   cat(paste0(col_start, text, col_end, "\n"))
}
catcyan <- function(text) {
   col_start <- "\033[0;0;36m"
   col_end <- "\033[0m"
   cat(paste0(col_start, text, col_end, "\n"))
}
catmagenta <- function(text) {
   col_start <- "\033[0;0;35m"
   col_end <- "\033[0m"
   cat(paste0(col_start, text, col_end, "\n"))
}



###### TREE PROCESSING

getLeafNames <- function(tree, taxa, suffix) {
   full_names <- tibble(
   						leaves=tree$tip.label,
   						name=tree$tip.label
   					)
   if (length(suffix) > 0) {
   	full_names <- tibble(
   							leaves=tree$tip.label,
   							name=sub(
   								paste0(suffix,"$"),
   								"",
   								tree$tip.label,
   								perl=TRUE
   								)
   							)
   }
   leaves <- inner_join(full_names, taxa, by='name')
   return(leaves)
}



taxaFromNames <- function(tree,taxa,separator,suffix) {

   if (separator == "|") {
      separator <- "\\|"
   }

   leaves <- tree$tip.label
   leaves <- grep(separator,leaves,value=TRUE)
   groupsFromNames <- tibble(
      name=sub(
				paste0(suffix,"$"),
				"",
				leaves,
				perl=TRUE
				),
      group=sub(
         paste0(separator,".*"),
         "",
         leaves,
         perl=TRUE
         )
      )

   if(length(taxa) > 0) {
      groupsFromNames_aj <- anti_join(groupsFromNames,taxa,by="name")
   }
   if(length(taxa) == 0) {
      taxa <- tibble(name=character(),group=character())
   }
   groups <- full_join(taxa,groupsFromNames,by=c("name","group"))

   return(groups)
}



calculate_x_limit <- function(tree) {

   t <- ggtree(tree)
   x_limit <- t$data %>% filter(isTip == TRUE) %>% select(x) %>% max()
   x_limit = x_limit + x_limit/3
   return(x_limit)
}

checkMono <- function(tree, leaves, node, group) {

   offspring_leaves <- getOffspringLabels(tree,node)
   offspring_groups <- leaves %>%
                           filter(leaves %in% offspring_leaves$label) %>%
                           select(group) %>%
                           unique()
   isGroup <- ifelse(length(offspring_leaves$label) >1, "TRUE", "FALSE")
   present <- ifelse(group %in% offspring_groups$group, "TRUE", "FALSE")
   monophyl <- ifelse(length(offspring_groups$group) == 1, "TRUE", "FALSE")

   decision <- ""
   ifelse(isGroup,
      ifelse(present,
         ifelse(monophyl,
            decision <- "monophyletic",
            decision <- "present"
         ),
         decision <- "absent"
      ),
      decision <- "monotypic"
   )
   return(decision)
}

checkChildrenMono <- function(tree, leaves, node, group, node_list) {

   children <- tree %>% child(node)

   for (i in 1:2) {
      result <- checkMono(tree, leaves, children[i], group)
      if (result == "monophyletic") {
         node_list <- c(node_list, children[i])
      }
      # if (result == "monotypic") {
      #    node_list <- c(node_list, children[i])
      # }
      if (result == "present") {
         node_list <- checkChildrenMono(tree, leaves, children[i], group, node_list)
      }
   }
   return(node_list)
}

getOffspringLabels <- function(tree, node) {
   offspring_labels <- tree %>%
                           as_tibble() %>%
                           offspring(node) %>%
                           select(label) %>%
                           filter( str_detect(label,"[A-Za-z]")) %>%
                           drop_na()
   return(offspring_labels)
}

monophyletic_subgroups <- function(tree, leaves, col_groups,threads) {

#   groupOTUs = list()
   groupMono = data.frame(group=1,node=1,size=1,col=1)

   groups <- leaves %>% select(group) %>% unique()
   groups <- groups$group %>% sort()

   groupOTUs <- sapply(groups, listGroupOTUs)

   if (threads == 1) {
      groupMono <- lapply(groups, groupAnalysis)
   } else {
      groupMono <- mclapply(groups, groupAnalysis, mc.cores=threads)
   }

   groupMono <- groupMono %>%
      reduce(full_join,by=c("group","node","size","col"))

   groupMono$node <- as.numeric(groupMono$node)
   groupMono$size <- as.numeric(groupMono$size)

   return_list <- list(groupMono, groupOTUs)
   return(return_list)
}


listGroupOTUs <- function(group) {
   g = group
   extract_leaves <- leaves %>% filter(group == g) %>% select(leaves)
   leafNames <- extract_leaves$leaves

   return(leafNames)
}


groupAnalysis <- function(group) {
    g <- group
    extract_leaves <- leaves %>% filter(group == g) %>% select(leaves)
    leafNames <- extract_leaves$leaves
    col_g <- col_groups %>% filter(group==g) %>% select(col)
    if (nrow(col_g) == 0) {
       col_g <- tibble(group=g, col="black")
    }

    groupMono <- data.frame(
      group=character(),
      node=numeric(),
      size=numeric(),
      col=character()
   )

    if (length(leafNames) > 1) {

       node <- MRCA(tree, leafNames)

       result <- checkMono(tree, leaves, node, g)
       if (result == "monophyletic") {
          lab <- getOffspringLabels(tree, node)
          n_group <- lab %>% nrow()

          groupMono <- rbind(groupMono, data.frame(group=g, node=node, size=n_group, col=col_g$col))
          catgreen(paste0("Collapsing ", g, "..."))
       }
       if (result == "present") {
          nodesMono <- vector()
          nodesMono <- checkChildrenMono(tree, leaves, node, g, nodesMono)
          n_nodesMono <- length(nodesMono)
          if (n_nodesMono >= 1) {
             for (n in 1:n_nodesMono) {
                lab <- getOffspringLabels(tree, nodesMono[n])
                n_group <- lab %>% nrow()
                # if (n_group == 0) {
                #    n_group <- 1
                # }
                groupMono <- rbind(groupMono, data.frame(group=g, node=nodesMono[n], size=n_group, col=col_g$col))

                catcyan(paste0("Collapsing ", g, ", ", n, "/", n_nodesMono,"..."))
             }
          }
          if (n_nodesMono == 0) {
             catmagenta(paste0("Not collapsing ", g, "..."))
          }
       }
    }
   return(groupMono)

}



###### TREE PLOTTING

collapse_tree <- function(plot, toCollapse, nodesTotal) {

   nodes <- length(toCollapse$node)

   for (i in 1:nodes) {
      size <- toCollapse[i,"size"]
      plot <- scaleClade(
         plot,
         node=toCollapse[i,"node"],
         scale=2/size
      )
   }
   for (i in 1:nodes) {
       g <- toCollapse[i,"group"]
       n <- toCollapse[i,"node"]
       s <- toCollapse[i,"size"]
       col <- toCollapse[i,"col"]
       label = paste0(g, " (", s, ")")

       d <- data.frame(node=n, name=g, col=col,lab=label)

       plot <- plot + geom_cladelab(
         data=d,
         mapping=aes(
   	       node=node,
   	       label=lab,
   	    ),
         fontsize=3,
         textcolour=col,
         barsize=0,
         barcolour="white"
       )
       plot <-  plot %>%
               collapse(
                     node=n,
                     'mixed',
                     fill=col,
                     col="black"
                  )
       plot <- plot + theme(legend.position="none")
   }

   ##### Suppresing warnings due to unclear warning message:
   ##### Removed x rows containing missing values (geom_point_g_gtree).
   oldw <- getOption("warn")
   options(warn=-1)
   return(plot)
   options(warn=oldw)
}

color_branches <- function(plot, groupedOTUs, color_vector_groups) {

   plot <- groupOTU(plot, groupedOTUs, 'group') +
        aes(color=group) +
        scale_color_manual(values=color_vector_groups) +
        theme(legend.position="none")
   return(plot)
}

calculate_thresholds <- function(plot, low, medium, high) {
   type <- ""

   max <- plot$data %>% filter(str_detect(label,"^[0-9]+$")) %>% select(label) %>% simplify() %>% as.numeric() %>% max()

   ifelse(max > 1,
      type <- "percent",
      type <- "fraction"
   )
   if (type == "fraction") {
      low <- low/100
      medium <- medium/100
      high <- high/100
   }

   thresholds=list("low"=low, "medium"=medium, "high"=high)
   return(thresholds)
}

draw_support_values <- function(plot, low_threshold=70, medium_threshold=95, high_threshold=100) {

   thresholds <- calculate_thresholds(plot, low_threshold, medium_threshold, high_threshold)

   plot$data <- plot$data %>% mutate(
      support.category =
         ifelse(
            isTip == TRUE,
            NA,
            ifelse(
               label == "",
               NA,
               ifelse(
                  suppressWarnings(as.numeric(label)) >= thresholds$high,
                  3,
                  ifelse(
                     suppressWarnings(as.numeric(label)) >= thresholds$medium,
                     2,
                     ifelse(
                        suppressWarnings(as.numeric(label)) >= thresholds$low,
                        1,
                        NA
                     )
                  )
               )
            )
         )
      )
   plot$data$support.category <- factor(plot$data$support.category, levels=c(1,2,3))

    eval(substitute(
       expr={
       plot <- plot +
            geom_nodepoint(
         	   aes(
                  x=branch,
                  subset =
         		   ! is.na(suppressWarnings(as.numeric(label))) &
         		   suppressWarnings(as.numeric(label)) >= thresholds$low,
                  fill = support.category
         	   ),
               color="black",
               shape=21,
         	   size=1.5
            ) +
            scale_fill_manual(
               name="Support values",
               labels=c(
                  paste0("[",thresholds$low,",",thresholds$medium,")"),
                  paste0("[",thresholds$medium,",",thresholds$high,")"),
                  paste0("≥",thresholds$high)
               ),
               values=c("white", "gray67", "black")
            ) +
            theme(legend.position="right")
    },
       env=list( thresholds=thresholds )
    ))
   return(plot)
}

print_support_values <- function(plot) {
   plot <- plot +
         geom_tiplab(size=1) +
         geom_nodelab(
   	      aes(
               x = branch,
               subset = ! is.na(suppressWarnings(as.numeric(label)))
   	      ),
        	   nudge_y=0.75,
   	      hjust=0.5,
        	   size=1
         )
   return(plot)
}

print_tiplabels <- function(plot, taxa) {
   plot$data$label <- gsub("_"," ",plot$data$label)
   plot$data$label <- gsub("GCA ","GCA_",plot$data$label)

   plot <- plot + geom_tiplab(
                        size=1.5
                     )
}

draw_root <- function(tree, plot, root)  {
   rootedge <- NULL

   if (length(tree$root.edge) == 1) {
      rootedge <- tree$root.edge
   }
   if (root && length(tree$root.edge) == 0) {
      rootedge <- plot$data %>%
                     select(branch.length) %>%
                     simplify() %>%
                     mean() %>%
                     `/`(2)
   }
   plot <- plot + geom_rootedge(rootedge = rootedge)

   return(plot)
}
