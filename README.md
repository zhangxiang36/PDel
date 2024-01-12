# Aanalyze outcome of gene-editing outcome by calling indel/SNP with NGS data

Scripts for running the analysis gene-editing outcome, such as Knock-out, base-editing, PE and deletion.

## Running the script

### Calling indel/SNP

To perform analysis of NGS data for indel/SNP, clone this repository and run
```bash

# data structure

└── NGS
    ├── D18-PDel-Syn
    │   ├── PY5704-H50-50-10.R1.fq.gz
    │   ├── PY5704-H50-51-5.R1.fq.gz
    │   └── Ref.txt
    ├── OsALS-PDel-Syn
    │   ├── PY5323-H50-22-14.R1.fq.gz
    │   ├── PY5323-H50-2-4.R1.fq.gz
    │   └── Ref.txt
    ├── OsGS3-PDel
    │   ├── PY5473-H50-01-1-15.R1.fq.gz
    │   ├── PY5473-H50-03-1-23.R1.fq.gz
    │   └── Ref.txt
    ├── OsGS3-PE3
    │   ├── PY4448-H50-20-47.R1.fq.gz
    │   ├── PY4448-H50-21-1-45.R1.fq.gz
    │   └── Ref.txt
    └── sample_annotation.xlsx


/path/to/SNP_indel.sh /path/to/NGS

```

Results will be shown in <sample_name>_all_sample_combined.txt (e.g., D18-PDel-Syn_all_sample_combined.txt). The bam files can be checked with IGV or Tablet software.

### Denovo assemble 

For deletions of hundreds of base pair, we use denovo  methods to check mutation outcome of gene-editing.

To perform denovo analysis, clone this repository and run:

```bash

# data structure

└── NGS
    └── OsGS3-504-PE3
         ├── PY5505-H50-02-51.R1.fq.gz
         ├── PY5505-H50-06-1.R1.fq.gz
         └── Ref.txt


/path/to/denovo.sh /path/to/NGS/OsGS3-504-PE3

```

Results will be stored in Trinity-GG.fasta and RSEM.isoforms.results for assmebled sequence and quality of assembled sequences, respectively.
