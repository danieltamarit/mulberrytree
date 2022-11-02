# mulberrytree

Produce publication-quality figures out of newick tree files and user-defined groups of leaves.

Usage
```
Rscript mulberrytree.R tree=<newick_tree> groups=<leaf_classification> colors=<group_colors>
```

Example run:
```
Rscript mulberrytree.R tree=example/a_treeSimple.nwk groups=example/b_taxa.tsv colors=example/c_colorGroups.tsv
```
