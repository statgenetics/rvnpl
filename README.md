# RV-NPL manual

## Installation

### Requirements

+ python: 2.7

+ numpy: >= 1.11.0

+ boost_python

  or

+ Anaconda: >= 2.3

The RV-NPL package is free and available on github.  Run the following commands to download and install the RV-NPL.

``` shell
git checkout https://github.com/percylinhai/rvnpl.git
cd rvnpl
python setup.py install 
```

If the program is installed correctly, you will see program options using the following command:

```shell
rvgdt --help
```



## Input Format

### Genotype File

The genotype file gives the genotype information of each subject per line. No header is needed in the genotype file. The first column is the subject id (which should be the same as the subject id in pedigree files), and the following columns is the number of minor allele on each variant sites (0/1/2 coding and -9 is missing), which is separated by a space or tab. An example of the genotype file is given below

```
11000.fa -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9
11000.mo -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9
11000.p1 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9
11000.s1 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9
11001.fa -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9
11001.mo -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9
11001.p1 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9
11002.fa -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9
11002.mo -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9
11002.p1 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9
```



### Pedigree File

The pedigree file is a white-space (space or tab) delimited file without header. The first six columns are mandatory:

+ Family ID     
+ Individual ID: must be unique within the family     
+ Paternal ID: 0 if not available 
+ Maternal ID: 0 if not available 
+ Sex:  1=male, 2=female   
+ Phenotype: 1=unaffected, 2=affected

An example pedigree file is given below:

```
11000 11000.fa 0 0 1 1
11000 11000.mo 0 0 2 1
11000 11000.p1 11000.fa 11000.mo 1 2
11000 11000.s1 11000.fa 11000.mo 2 1
11001 11001.fa 0 0 1 1
11001 11001.mo 0 0 2 1
11001 11001.p1 11001.fa 11001.mo 1 2
11002 11002.fa 0 0 1 1
11002 11002.mo 0 0 2 1
11002 11002.p1 11002.fa 11002.mo 2 2
```

### Optional Files

#### Selected variant file

In case you want to weight each variant differently. The weights can be given in a single-column file (no header), in which each line is the weight for the corresponding variant site (the order should be the same as the order of variant sites in genotype file).

## Options

### Options for generating CHP markers

```
TBD
```

Example commands are shown below:

```shell
rvnpl collapse --fam 100extend_01.ped --vcf $datapath/$g/rep1.vcf.gz -f MERLIN --output ./rep1 --freq EVSMAF -c 0.01 --rvhaplo --include_vars $g.txt 

rvnpl collapse --fam 100extend_01.ped --vcf $datapath/$g/rep1.vcf.gz -f MERLIN --output ./rep1 --freq_by_fam Freq_by_Fam.txt -c 0.01 --rvhaplo --include_vars $g.txt
```

### Options for npl analysis

```
optional arguments:
  -h, --help            show this help message and exit

Options for doing analysis on CHP markers(default) or SNV markers:
  --snv                 Calculate on SNV markers

Input/Output options:
  --path PATH           Path for input pedigree information.
  --output PATH         Path for output files
  --n_jobs N            number of multiprocess

Options for calculating p-values:
  --exact               get the exact distribution of Z-score and calculate p
                        value from it
  --cut FLOAT, -c FLOAT
                        threshold for adaptive permutations
  --rep N               times of permutations
  --fam_rep N           times of permutations for each family
  --force               keep permutation times unchanged
  --perfect_max N       maximum for inheritance vector iterations
  --info_only           include only informative families
  --perfect             use perfect data approximation in calculating Z-score
  --kc                  use Kong&Cox(1997) extension for analytical p-values
  --sall                enable calculation of NPL-all
  --rvibd               calculate IBD for RV only
```

Example commands are shown below:

```shell
rvnpl npl --path CHP_func_SEQL/${GENE} --output CHP_func_SEQL/${GENE} --exact --info_only --perfect --sall --rvibd --n_jobs 8 -c 0.001 --rep 2000000

```

The output is given in the ${proj}.rvgdt_output file. 

# Questions

If you have any further questions, please fell free to create a issue ticket in github. 
