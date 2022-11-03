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
mulberrytree -h
mulberrytree -t <newick_tree> -g <leaf_classification> -c <group_colors>
mulberrytree -t <newick_tree> -l -s "<separator>"
```

Example runs:
```
mulberrytree -t example/a_treeSimple.nwk -g example/b_taxa.tsv -c example/c_colorGroups.tsv -x "_[0-9]+G"
mulberrytree -t example/a_treeSimple2.nwk -l
mulberrytree -t example/a_treeSimple3.nwk -l -s "___"
```
