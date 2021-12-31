# Sankoff is a fast and parallel implementation for ancestral state reconstruction based on Sankoff's algorithm
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/sankoff/README.html)

# Installation
**Conda installation (Linux and Mac):**
(Generated based on README from [CRISPRme](https://github.com/samuelecancellieri/CRISPRme))
- Open a terminal window
- Paste this command into the terminal:
    ```
    curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh --output Miniconda3-latest-Linux-x86_64.sh
    ```
- If the file is correctly downloaded you now need to execute it to complete the installation, so paste this command into the terminal:
    - Linux
    ```
    bash Miniconda3-latest-Linux-x86_64.sh
    ```
- Press ENTER when requested and yes when an answer is requested, in this way you allow conda to set all the directories in your HOME path for an easy use
- After the complete installation you will receive this message “Thank you for installing Miniconda3!” to certify the correct installation.
- Now you need to close the terminal window you are using and open a new one, to allow the system to start conda.
- In the new terminal window you should see something like this:
    ```
    (base) user@computer:~$
    ```
    If you read the "(base)" like this, conda is loaded correctly and you can start using it.
- Now, you can install sankoff by typing the command:
    ```
    conda install -c bioconda sankoff
    ```
- To test your installation, type the command:
    ```
    sankoff --help
    ```
- After the execution of the command you should see sankoff options.

Now the software is installed and ready to be used.


# Usage
An easy usage could be as follows
```
sankoff --tree [TREE-FILE] --aln [INPUT-SEQUENCE-FILE] --cost [COST-MATRIX-FILE] --out-as [ANCESTRAL-STATE-OUTPUT-FILE]
```
Sankoff with the above parameters reads the tree from `TREE-FILE` and leaf node sequences from `INPUT-SEQUENCE-FILE` as well as cost matrix from file `COST-MATRIX-FILE`, executes the sankoff algorithm and prints sequence for leaf nodes and internal nodes to file `ANCESTRAL-STATE-OUTPUT-FILE`.
If internal nodes are unnamed, option `--ilabel [PREFIX]` could be used which assigns label to internal nodes with prefix `PREFIX`. Then the tree with new names could be saved with option `--out-tree [OUT-TREE]` in file `OUT-TREE`.

## Example
  - Download leaf node sequences from [aln.txt](https://raw.githubusercontent.com/hzi-bifo/sankoff/main/tests/data/aln1.txt) and save as aln1.txt
  - Download cost matrix from [cost-1.txt](https://raw.githubusercontent.com/hzi-bifo/sankoff/main/tests/data/cost-1.txt) and save as cost-1.txt
  - Download input phylogeny from [tree1.txt](https://raw.githubusercontent.com/hzi-bifo/sankoff/main/tests/data/tree1.txt) and save as cost-1.txt
  - Run 
```
sankoff --tree tree1.txt --aln aln.txt --cost cost-1.txt --ilabel inode --out-as aln-out.txt --out-tree tree-out.txt
```
  - Internal nodes of the input tree is now labled and the tree is saved in `tree-out.txt`.
  - Sequence infered for internal nodes as well as sequence for leaf nodes are saved in `aln-out.txt`.

# Options
  * `--nthread [THREAD-COUNT]`: Use `THREAD-COUNT` threads for execution.
  * `--tree [TREE-FILE]`: Input phylogeny in the plain tree format.
  * `--nexus [TREE-FILE]`: Input phylogeny in the nexus format.
  * `--aln [INPUT-SEQUENCE-FILE]`: Input sequences for leaf nodes.
  * `--cost [COST-MATRIX-FILE]`: File containing the cost matrix. First row and first column represent characters in the sequence files. Values are integer.
  * `--cost-identity-aa [identical] [non-identical] [X-AA] [X-X] [X-GAP] [AA-GAP] [GAP-GAP]`: Instead of cost matrix, this option may be used. This option builds a cost matrix for amino acids with two special characters: gap (`-`) and unknown (`X`). The matrix is filled with the parameters given to this option.
  * `--cost-identity-dna [identical] [non-identical] [X-NA] [X-X] [X-GAP] [NA-GAP] [GAP-GAP]`: Instead of cost matrix, this option may be used. This option builds a cost matrix for DNA with two special characters: gap (`-`) and unknown (`X`). The matrix is filled with the parameters given to this option.
  * `--out-as [ANCESTRAL-STATE-OUTPUT-FILE]`: Output file containing all the sequences for internal nodes as well as leaf nodes.
  * `--ilabel [PREFIX]`: If present, internal nodes are labeled with `[PREFIX]#NUM` where `#NUM` will be different values.
  * `--out-tree [OUT-TREE]`: The phylogeny, potentially after labeling, will be saved in this file.
  * `--induce_tree_over_samples`: Remove samples from the tree which are not present in the sequence file and rebuild the tree.
  
# Contributors
  * Mohammad-Hadi Foroughmand-Araabi
  * Sama Goliaei

# License
This project is under GPL 3.0 license.
