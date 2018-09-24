# Misc scripts

## Match and Symlink script
This script is used to match specific files in a directory and create symbolic links of the matched files to another directory. Handy if you have to symlink specific files in a directory with many other files.

Usage:

**Rscript match_and_symlink.R file_list file_dir output_dir pattern 
recursive_setting**


- file_list: A list of file names to match, tab separated file without 
a header.

- file_dir: The directory of the files in file_list.

- output_dir: The directory where the symlinked files will be placed.

- pattern: The pattern to look for, f. ex. "fastq.gz" or "zip".

- recursive_setting: TRUE/FALSE, specifies if the script should look for 
files in sub-folders.

## Fastqc analysis script
R script for analysis of fastqc reports. The script relies on the .zip
files created by fastqc, and use the fastqcr package
(https://CRAN.R-project.org/package=fastqcr) for importing the data (no
need to unpack the .zip files).

The script creates a few informative plots based on the data from fastqc
and saves it in a new folder called "fastqc_results_TODAYS_DATE" in the
output_dir location.

Usage:

**Rscript fastqc_analysis.R zipfiles_location output_dir**
