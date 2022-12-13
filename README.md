# PostKraken

A set of tools to process the output from Kraken2

## About:

### Motivation:

Hello and welcome I am a PhD student in the BBSRC LIDo program; while doing my rotations, I was shocked to find out how few tools there are available to process the output of Kraken2, and the tools that did exist didn't seem to work very well for my use case.

Thus I endeavoured to write my own scripts to do the tasks I needed and decided to share them for others who might also find them useful ...(dramatic music)... and PostKraken was born!

### Plans moving forward:

I intend this to very much be an active project under development. My current plan is to add tools and features as I need them for my work. But if there are any tools or features you would like to see added, please make yourself known, and I'll do my best to accommodate you.

## Dependencies:

In order to run any of the scripts here you'll need to have the following installed on your PATH:

-   Python (version 3.2 or higher)

-   R

    -   data.table

    -   stringi

    -   argparse

## Running Scripts:

All scripts include here are ready to go no install required, all you need to do is to download or copy the scripts into a local file and run the following commands and your ready to rock and roll

```{bash}
chmod +x ${pathToScript}/Script.R
${pathToScript}/Script.R -h
```

## filtering_kraken_output.R:

This script is used to filter through the output of kraken2 and extract/exclude reads identified as a specific taxa(s)

### Usage and arguments:

```{bash}
usage: filtering_kraken_output.R [-h] [{PE,SE,C}] -k KRAKEN_FILE [-k2 KRAKEN_FILE_2]
                                  -o OUTPUT_FILENAME [-o2 OUTPUT_FILENAME_2] 
                                  -m [{filter,exclude}] -t TAXA_ID 
                                  [--include_lower_taxa] [{True,False}]]
                                  [-r KRAKEN_REPORT_FILENAME]


positional arguments:
  {PE,SE,C}             A flag to specify the input type: PE = paried end, SE = single end, 
                        C = assembly file containing contigs

optional arguments:
  -h, --help                          show this help message and exit
  
  -k KRAKEN_FILE, -k1 KRAKEN_FILE, --kraken_output KRAKEN_FILE
                                      The file name of the kraken output file you want to
                                      filter
                                      
  -k2 KRAKEN_FILE_2, --kraken_output_2 KRAKEN_FILE_2
                                      The file name of the kraken output file for your
                                      reverse reads
                                      
  -o OUTPUT_FILENAME, -o1 OUTPUT_FILENAME, --output_filename OUTPUT_FILENAME
                                      The file name you want to give the filtered reads
                                      
  -o2 OUTPUT_FILENAME_2, --output_filename_2 OUTPUT_FILENAME_2
                                      The file name you want to give the filtered reverse
                                      reads
                                      
  -m [{filter,exclude}], --method [{filter,exclude}]
                                      filter finds reads that match specified taxa IDs.
                                      exclude finds reads that do not match the specified
                                      taxa ID
                                      
  -t TAXA_ID, --taxa TAXA_ID
                                      The taxa ID you want to filter/exclude by
                                      
  --include_lower_taxa [{True,False}]
                                      whether to include lower taxa or not. Defualt is F if
                                      the T value is passed you will need to include the -r
                                      flag for it to work
                                      
  -r KRAKEN_REPORT_FILENAME, --kraken_report KRAKEN_REPORT_FILENAME
                                      The report file from your kraken run
```

### Examples:

In this section, you will find a series of examples of how to use the flags outlined to achieve desired outcomes. In all of these examples, we will use the taxonomic ID associated with the genus Yersinia (taxaID: 629)

#### 1) Filtering paired end data for a supplied taxonomic ID:

```{bash}
./filtering_kraken_output.R PE -k F_Reads_classified.fasta -k2 R_Reads_classified.fasta -o F_Reads_filtered.fasta -o2 F_Reads_filtered.fasta -m filter -t 629 
```

#### 2) Excluding a supplied taxonomic ID from paired end data:

```{bash}
./filtering_kraken_output.R PE -k F_Reads_classified.fasta -k2 R_Reads_classified.fasta -o F_Reads_filtered.fasta -o2 F_Reads_filtered.fasta -m exclude -t 629 
```

#### 3) Filtering an assembly for contigs assigned to a supplied taxonomic ID and those assigned at a more specific level:

```{bash}
./filtering_kraken_output.R C -k assembly_classified.fasta -o assembly_filtered.fasta -m filter -t 629 --include_lower_taxa T -r assembly_report
```

#### 4) Excluding contigs assigned to a supplied taxonomic ID and those assigned at a more specific level from an assembly:

```{bash}
./filtering_kraken_output.R C -k assembly_classified.fasta -o assembly_filtered.fasta -m exclude -t 629 --include_lower_taxa T -r assembly_report
```
