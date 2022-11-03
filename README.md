# mulberrytree

Produce publication-quality figures out of newick tree files and user-defined groups of leaves.


### Install dependencies
```
# Clone repository
git clone https://github.com/danieltamarit/mulberrytree.git

# Install basic dependencies
conda env create -f envs/conda_env_base.yml

# Install ggtree through Bioconductor
R
> if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

> BiocManager::install("ggtree")
> quit()
```


### Usage
```
Rscript mulberrytree.R tree=<newick_tree> groups=<leaf_classification> colors=<group_colors>
```

Example run:
```
Rscript mulberrytree.R tree=example/a_treeSimple.nwk groups=example/b_taxa.tsv colors=example/c_colorGroups.tsv
```
