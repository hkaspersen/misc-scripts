# Misc scripts

## Match and Symlink script
This script is used to match specific files in a directory and create symbolic links of the matched files to another directory. Handy if you have to symlink specific files in a directory with many other files.

Usage:

```
Rscript match_and_symlink.R file_list file_dir output_dir pattern 
recursive_setting
```

- **file_list**: A list of file names to match, tab separated file without a header.

- **file_dir**: The directory of the files in file_list.

- **output_dir**: The directory where the symlinked files will be placed.

- **pattern**: The pattern to look for, f. ex. "fastq.gz" or "zip".

- **recursive_setting**: TRUE/FALSE, specifies if the script should look 
for files in sub-folders.

## Fastqc analysis script
R script for analysis of fastqc reports. The script relies on the .zip
files created by fastqc, and use the fastqcr package
(https://CRAN.R-project.org/package=fastqcr) for importing the data (no
need to unpack the .zip files).

The script creates a few informative plots based on the data from fastqc
and saves it in a new folder called "fastqc_results_TODAYS_DATE" in the
output_dir location.

Usage:

```
Rscript fastqc_analysis.R zipfiles_location output_dir
```

## Mash screen analysis script

Script for visualizing and identifying contaminants in read files, based on results from Mash screen. 
The script creates three files in the specified output directory: One plot and two text files. 

Usage:

```
Rscript mash_screen.R report_dir file_pattern filter_value organism output_dir
```

- **report_dir**: Full path to the directory of the mash screen reports

- **file_pattern**: Suffix of mash screen reports (ex. .out, _out, etc.)

- **filter_value**: A value set for filtering out results below given 
value, based on shared hashes. Value between 0-1.

- **organism**: name of the organism of interest, for example 
"escherichia coli". Include quotation marks if two words are used. This organism is presumed to be the "wanted" organism in the files.

- **output_dir**: The directory where the output files are deposited.


Output files:

- **contaminated_samples_report.txt**: A tab separated text file 
containing the samples with possible contamination. These are filtered 
out based on the presence of Identity values above 0.95 other than the 
organism of interest.

- **full_mash_report.txt**: A tab separated text file containing the 
mash screen results where shared hashes were 100/1000 or higher. 

- **mash_plot.svg**: A dotplot that visualizes the results from the mash 
screen analysis. A cut-off value for the identity value is presented as 
a horizontal line, where all dots above this line should optimally be of 
the organism of interest. If not, significant contamination is likely.

## Extract gene sequences script
Script for extracting nucleotide or amino acid sequences from multiple 
.ffn or .faa files created by prokka. Outputs a fasta file with the gene 
sequences.

Usage:
```
Rscript extract_gene.R fasta_file_location gene_name nt_or_aa output_dir
```
- **fasta_file_location**: The folder where the prokka results are located
- **gene_name**: Name of the gene, as listed in the .faa or .ffn fasta headers
- **nt_or_aa**: Either "aa" for amino acid sequences or "nt" for nucleotide sequences
- **output_dir**: Location of output file


## Create manifest files for ENA submission
Script for creating separate manifest files for each sample, used in ENA submission of fastq files.
Input is a list over all fastq file names in the submission.

Usage:
```
Rscript create_manifest_file.R list_file output_dir study name instrument library_name library_source library_selection description
```

See https://ena-docs.readthedocs.io/en/latest/reads/webin-cli.html for details on values.