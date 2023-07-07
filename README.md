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

-   R

    -   data.table

    -   stringi

    -   optparse

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
Usage: filtering_kraken_output.R --i [PE,SE,C] --k KRAKEN_FILE [--k2 KRAKEN_FILE_2] 
                                 --o OUTPUT_FILENAME [--o2 OUTPUT_FILENAME_2]
				                         --m [F,E] --t TAXA_ID [--include_lower_taxa [T,F]] 
				                         [--r KRAKEN_REPORT_FILENAME]

Description: Extracting reads identified as a specific taxa(s)

Options:
	--input_type=[PE, SE, C]
		A flag to specify the input type: PE = paired end, SE = single end, C = assembly file     
		containing contigs, this is a required flag

	--i=[PE, SE, C]
		A flag to specify the input type: PE = paired end, SE = single end, C = assembly file 
		containing contigs, this is a required flag

	--kraken_output=[KRAKEN2 OUTPUT FILE]
		The file name of the kraken output file you want to filter, this is a required flag

	--k=[KRAKEN2 OUTPUT FILE]
		The file name of the kraken output file you want to filter, this is a required flag

	--kraken_output_2=[KRAKEN2 OUTPUT FILE]
		The file name of the kraken output file for your reverse reads

	--k2=[KRAKEN2 OUTPUT FILE]
		The file name of the kraken output file for your reverse reads

	--output_filename=[NAME OF OUTPUT FILE]
		The file name you want to give the filtered reads

	--o=[NAME OF OUTPUT FILE]
		The file name you want to give the filtered reads

	--output_filename_2=[NAME OF OUTPUT FILE]
		The file name you want to give the filtered reverse reads

	--o2=[NAME OF OUTPUT FILE]
		The file name you want to give the filtered reverse reads

	--method=[F, E]
		F = filter, E = exclude. Filter finds reads that match specified taxa IDs.
		Exclude finds reads that do not match the specified taxa ID, this is a required flag

	--m=[F, E]
		F = filter, E = exclude. Filter finds reads that match specified taxa IDs.
		Exclude finds reads that do not match the specified taxa ID, this is a required flag

	--taxa=[TAXA ID]
		The taxa ID you want to filter/exclude by, this is a required flag

	--t=[TAXA ID]
		The taxa ID you want to filter/exclude by, this is a required flag

	--include_lower_taxa=[F, T]
		whether to include lower taxa or not: F = do not include lower taxa, T = include lower 
		taxa.Defualt is F if the T value is passed you will need to include the --r flag for it
		to work

	--r=[KRAKEN2 REPORT FILE]
		The report file from your kraken run

	--kraken_report=[KRAKEN2 REPORT FILE]
		The report file from your kraken run

	-h, --help
		Show this help message and exit
```

### Examples:

In this section, you will find a series of examples of how to use the flags outlined to achieve desired outcomes. In all of these examples, we will use the taxonomic ID associated with the genus Yersinia (taxaID: 629)

#### 1) Filtering paired end data for a supplied taxonomic ID:

```{bash}
./filtering_kraken_output.R --i PE --k F_Reads_classified.fasta --k2 R_Reads_classified.fasta --o F_Reads_filtered.fasta --o2 R_Reads_filtered.fasta --m F --t 629 
```

#### 2) Excluding a supplied taxonomic ID from paired end data:

```{bash}
./filtering_kraken_output.R --i PE --k F_Reads_classified.fasta --k2 R_Reads_classified.fasta --o F_Reads_filtered.fasta --o2 R_Reads_filtered.fasta --m E --t 629 
```

#### 3) Filtering an assembly for contigs assigned to a supplied taxonomic ID and those assigned at a more specific level:

```{bash}
./filtering_kraken_output.R --i C --k assembly_classified.fasta --o assembly_filtered.fasta --m F --t 629 --include_lower_taxa T --r assembly_report
```

#### 4) Excluding contigs assigned to a supplied taxonomic ID and those assigned at a more specific level from an assembly:

```{bash}
./filtering_kraken_output.R --i C --k assembly_classified.fasta --o assembly_filtered.fasta --m E --t 629 --include_lower_taxa T --r assembly_report
```
