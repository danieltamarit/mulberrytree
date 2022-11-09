# mulberrytree

Produce publication-quality figures out of newick tree files and user-defined groups of leaves.


### Install dependencies
```
# Clone repository
git clone https://github.com/danieltamarit/mulberrytree.git
cd mulberrytree

# Install basic dependencies
conda env create -f envs/conda_env_base.yml
conda activate mulberry

```


### Options
```
-h     Print this help
-t     Input tree in Newick format (required)
-g     Group information in tsv format
         (Column 1: Leaf name; Column 2: Group name)
-c     Color information in tsv format
         (Column 1: Group name; Column 2: R-readable color)
-m     Midpoint root
-o     Prefix for output files
-l     Interpret group name from leaf name
-s     Separator for group interpretation from leaf (default: \"|\")
-x     Text or regular expression to be ignored as leaf name suffix
-T     Number of threads used for tree processing (default: 1)
         Note: optimal speed often reached with 1-2 threads
```

### Basic usage
```
mulberrytree -h
mulberrytree -t <newick_tree> -g <leaf_classification> -c <group_colors>
mulberrytree -t <newick_tree> -l -s "<separator>"
```

### Example runs:
```
mulberrytree -t example/a_treeSimple.nwk -g example/b_taxa.tsv -c example/c_colorGroups.tsv -x "_[0-9]+G"
mulberrytree -t example/a_treeSimple2.nwk -l
mulberrytree -t example/a_treeSimple3.nwk -l -s "___"
```



The name "mulberrytree" is a homage to figtree, a software I have used too many hours of my life using, and for which I often hoped I could have a less manual equivalent.
