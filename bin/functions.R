#/usr/bin/Rscript


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



###### TREE PROCESSING

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

monophyletic_subgroups <- function(tree, leaves, colortable) {

   groupOTUs = list()
   groupMono = data.frame(group=1,node=1,size=1,col=1)

   groups <- leaves %>% select(group) %>% unique()
   groups <- groups$group %>% sort()

   for (i in 1:length(groups)) {
       g = groups[i]
       extract_leaves <- leaves %>% filter(group == g) %>% select(leaves)
       leafNames <- extract_leaves$leaves
       col_g <- colortable %>% filter(group==g) %>% select(col)

       groupOTUs[[i]] <- leafNames
       names(groupOTUs)[i] <- g


       if (length(leafNames) > 1) {

          node <- MRCA(tree, leafNames)

          result <- checkMono(tree, leaves, node, g)
          if (result == "monophyletic") {
             lab <- getOffspringLabels(tree, node)
             n_group <- lab %>% nrow()

             groupMono <- rbind(groupMono, c(g, node, n_group, col_g$col))
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
                   groupMono <- rbind(groupMono, c(g, nodesMono[n], n_group, col_g$col))

                   catcyan(paste0("Collapsing ", g, ", ", n, "/", n_nodesMono,"..."))
                }
             }
             if (n_nodesMono == 0) {
                catcyan(paste0("Not collapsing ", g, "..."))
             }
          }
       }
   }
   groupMono <- groupMono[-1,]
   groupMono$node <- as.numeric(groupMono$node)
   groupMono$size <- as.numeric(groupMono$size)

   return_list <- list(groupMono, groupOTUs)
   return(return_list)
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
