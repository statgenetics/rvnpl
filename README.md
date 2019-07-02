# RV-NPL manual

## Installation

### Requirements

+ python: 2.7

+ Anaconda: >= 2.3

+ gcc-5, g++-5

+ boost_python

The RV-NPL package is free and available on github.  Run the following commands to download and install the RV-NPL.

``` shell
git clone https://github.com/percylinhai/rvnpl.git
cd rvnpl
python setup.py install 
```

If the program is installed correctly, you will see program options using the following command:

```shell
rvnpl --help
```



## Input Format

### Input files for generating CHP markers

#### Pedigree File

The pedigree file (PED file) is a white-space (space or tab) delimited file without header. The first six columns are mandatory:

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
#### Zipped and tabixed VCF file
The VCF file should contain variants for individuals corresponding to the PED file.

```
bgzip ./rep1.vcf
tabix -p vcf ./rep1.vcf.gz
```





### Optional Files

#### Selected variant file

If you want to analyze only subset of variants into CHP markers, you can provide a selected variant file with each row representing one variant (chromosome and position). For example:

```
19 58858740
19 58858782
19 58858787
```

#### Frequency by family file

When analyzing families from different ethnic populations, you need to provide this file to indicate the frequency column for each family. The frequency column should be the corresponding value under "INFO" in the VCF file. For example:

```
11000 gnomAD_NFE
11001 gnomAD_NFE
11002 gnomAD_AMR
```

## Options

### Options for generating CHP markers

```
optional arguments:
  -h, --help            show this help message and exit

Collapsed haplotype pattern method arguments:
  -b FILE, --blueprint FILE
                        Blueprint file that defines regional marker (format:
                        "chr startpos endpos name avg.distance male.distance
                        female.distance").
  --single-markers      Use single variant markers. This switch will overwrite
                        "--bin" and "--blueprint" arguments.

Input / output options:
  --fam FILE            Input pedigree and phenotype information in FAM
                        format.
  --vcf FILE            Input VCF file, bgzipped.
  --freq INFO           Info field name for allele frequency in VCF file.
  --freq_by_fam FILE    Per family info field name for allele frequency in VCF
                        file.
  --mle                 Estimate allele frequency from sample
  --rvhaplo             Only using rare variants for haplotyping
  -c P, --maf-cutoff P  MAF cutoff to define "common" variants to be excluded
                        from analyses.
  --include_vars FILE   Variants to be included in CHP construction
  --chrom-prefix STRING
                        Prefix to chromosome name in VCF file if applicable,
                        e.g. "chr".
  -o Name, --output Name
                        Output name prefix.

Runtime arguments:
  -j N, --jobs N        Number of CPUs to use.
  -q, --quiet           Disable the display of runtime MESSAGE.
```

Example commands are shown below:

```shell
cd example
rvnpl collapse --fam 100extend_01.ped --vcf A1BG/rep1.vcf.gz --output ./rep1 --freq EVSMAF -c 0.01 --rvhaplo --include_vars A1BG.txt 
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
rvnpl npl --path ./rep1 --output ./rep1 --exact --info_only --perfect --sall --rvibd --n_jobs 8 -c 0.001 --rep 2000000

```

The output is located in the given folder.

# Questions

If you have any further questions, please fell free to create a issue ticket in github. 
