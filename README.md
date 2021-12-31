# sankoff
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/sankoff/README.html)

Sankoff is a fast and parallel implementation for ancestral state reconstruction based on Sankoff's algorithm.

# Installation
**Conda installation (Linux and Mac):**
(Generated based on README from [CRISPRme(https://github.com/samuelecancellieri/CRISPRme)])
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
* TODO: download files and how to execute.

# Options
  * `--nthread [THREAD-COUNT]`: use `THREAD-COUNT` threads for execution.
  * `--tree [TREE-FILE]`
  * `--aln [INPUT-SEQUENCE-FILE]`
  * `--cost [COST-MATRIX-FILE]`
  * `--out-as [ANCESTRAL-STATE-OUTPUT-FILE]`
  * `--ilabel [PREFIX]`
  * `--out-tree [OUT-TREE]`
  
