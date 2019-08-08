# KRISPR -- CRISPR offtarget and efficiency prediction with *k*-mers

## Installation

First of all you have to clone the repository to your desired location. Then there are two ways to install all needed dependencies. But keep in mind, this tool only works on Linux/UNIX Systems.

### Manual Installation
Requirements are:
* Python 3 with:
    * numpy
    * kPal
    * pandas
    * sciPy
    * biopython
* jellyfish (can be downloaded from [here](https://github.com/gmarcais/Jellyfish))
* R (Version 3.2 or above)

All Python dependencies can be installed via pip:
```
pip install numpy kPal pandas scipy biopython
```

Make shure that all tools are available through the shell (appended to your `$PATH` Variable)

### Installation via Conda Environment (recommended)
The easier and more comfortable way is to use a Conda Environment. Hereby the only requirement is a working installation of [Anaconda](https://anaconda.org/) or [miniconda](https://conda.io/miniconda.html).

Just do:
```
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda create -n krispr python=3.6 numpy pandas scipy kPal biopython r
conda install -n krispr -c bioconda jellyfish
```
When you want to use the tool, you have to activate the environment by typing `source activate krispr`. When finished deactivate the environment by `source deactivate` to get back to your normal shell.

## Usage
`python krispr.py -h` gives you all the information you need.

## Examples

First you have to generate the *k*-mer Index with jellyfish by doing:

```
jellyfish count -m 21 -t 10 -o KINDEX.jf DATA.fastq
```

This builds the index with *k*-mer length 21, using 10 threads.

### Single target
Example call for a single target, coverage is 10 (`-c`), mismatches are 0 (`-m`):

```
python3.5 krispr.py single -q TCTGAGCCTCGGACCGACGGGGGG -j KINDEX.jf -m 0 -c 10
```

The output can look something like this:

```
This is KRISPR -- CRISPR offtarget and efficiency prediction with k-mers

Disclaimer: These values are only predictions and can not be taken as ground truth.

target:                    TCTGAGCCTCGGACCGACGGGGGG
mismatches:                0
off-target score :         5
GC all:                    0.75
GC distal:                 0.5556
GC proximal:               0.8571
entopy:                    1.2368
complexity:                0.2188
starts with A or G:        no
X on PAMX ends not with G: no

estimated efficiency: 0.7387
```

Please note that the input MUST be in the form: target + PAM + 1bp

### Multiple targets

Example call for the analysis of a FASTA file and the search for all potential targets, coverage is 10, mismatches are 0:

```
python3.5 krispr.py multi -q some_sequence.fasta -j KINDEX.jf -m 0 -c 10
```

The output then looks something like:

```
This is KRISPR -- CRISPR offtarget and efficiency prediction with k-mers

Disclaimer: These values are only predictions and can not be taken as ground truth.

seqid           target                          score   GC_all  GC_dist GC_prox ent complexity start_AG PAMX_G  e_eff
AY750996.1      TTCTCCCCCCCAATCCGCCCTGGG        0       0.7083  0.6667  0.7143  1.1646  0.2188  no      no      1.0
AY750996.1      CGATTTCCTCCGCGCCGTTCCGGT        0       0.6667  0.4444  0.7857  1.2031  0.1875  no      yes     0
AY750996.1      TCCTCCGCGCCGTTCCGGTCCGGC        0       0.7917  0.7778  0.7857  1.0327  0.1562  no      yes     0.5604
AY750996.1      CGCGCCGTTCCGGTCCGGCGAGGC        0       0.8333  0.7778  0.8571  1.1219  0.1406  no      yes     0.0986
AY750996.1      TCCGGTCCGGCGAGGCGACCCGGC        9       0.8333  0.7778  0.8571  1.1437  0.1094  no      yes     0
AY750996.1      GTCCGGCGAGGCGACCCGGCCGGC        10      0.875   0.7778  0.9286  1.0618  0.125   yes     yes     0
AY750996.1      GAGGCGACCCGGCCGGCCCCAGGG        10      0.875   0.7778  0.9286  0.9823  0.125   yes     no      0
AY750996.1      ATTTCCGCGCTGCTGCTGCTTGGT        3       0.5833  0.5556  0.5714  1.219   0.1562  yes     yes     0.5455
AY750996.1      CTGCTGCTGCTTGGTTCATCCGGT        0       0.5833  0.6667  0.5     1.219   0.1719  no      yes     0.0556
AY750996.1      TGCTCGCTCGCTGTTTCCTTCGGA        0       0.5833  0.6667  0.5     1.213   0.1562  no      yes     0.097
...
...
```

### Making a new model on your own data

The model shipped with the tool includes 32 targets from barley, wheat, tomato and tobacco.
With the following call it is possible to make a new model on your own data provided in the `targets.csv` file:

```
python3.5 krispr.py new-model -j KINDEX.jf -e targets.csv -c 10
```

The table file should be a `.csv` file and have a format like this:

```
target,efficiency
TTCTCCCCCCCAATCCGCCCTGGG,0.37
CGATTTCCTCCGCGCCGTTCCGGT,0.2
TCCTCCGCGCCGTTCCGGTCCGGC,0.0
CGCGCCGTTCCGGTCCGGCGAGGC,0.78
TCCGGTCCGGCGAGGCGACCCGGC,0.3
...
...
```

The tool reates a backup of the old file and then analyses all targets for their properties and writes a new data file.
Please consider that there should be at least 15-20 targets present to get sufficient predictions.


### Further explanations

In addition one can specify a parameter `-l` which indicates the CRISPR target length excluding PAM (defaults to 20). The minimum probe length in multi mode then has to be the `target length + PAM + 1` which in this example means 24. In single mode it has to be exactly `target length + PAM + 1`.
If one of the sequences in the FASTA in multi mode is to short krispr produces `--` in the output columns for this seqid. If the length fits, but there are no targets found it produces `-` in the specific output line.

Explanation of the columns:
1. seqid - name of the sequence from the FASTA
2. target - sequence of the found target(s) in the form: target + PAM + 1
3. score - off-taget score - normalized sum of the *k*-mer fequences for the target; the smaller this is, the better the target is, because there are less possible off-targets; good values are usually below 10 for 0 or 1 mismatch and below 50-150 for 2 or 3 mismatches
4. GC_all - overall GC content, should be between 0.5 to 0.8
5. GC_dist - distal (first 13 nucleotides) GC content
6. GC_prox - proximal (second half of the taget) GC content, distal and proximal GC content should not deviate significantly
7. ent - Shannon entopie of the target, should be around 0.6 +- 0.1
8. complexity - triplet diversity/complexity, how many different triplets are present in the target sequence regards to its length, should be between 1.2 to 1.3 +- 0.1
9. start_AG - does the target start with A or G
10. PAMX_G - is the nucleotide after the PAM not a G
11. e_eff - estimated efficiency of the target, calculated with a regression model taken all GC contents, entropie and complexity in to account, it ranges from 0 to 1 were 0 is the worst and 1 is the best. This value translates roughly to mutation frequencies.

To evaluate a sequence the most important values are the off-target score and the estimated efficiency.

### Info
KRISPR is part of the Kmasker tool (see [here](https://github.com/tschmutzer/kmasker)), so consider using this as well.
To cite this: ...
