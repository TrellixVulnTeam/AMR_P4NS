# AMR

The software you can find in this directory was used to simulate read liybraries and perform metagenomic analysis on them in order to identify antimicrobial resistance (AMR). In order to have this software should be enought to clone this repository in an envoirement where Anaconda is installed. The CAMISIM repository is a clone of the original camisim repo (https://github.com/CAMI-challenge/CAMISIM) that was modified in order to adchieve our golas.
To explain the functioning of the set of tools that was used is usefull to devide the process in three:

- Read library simulation
- Allignemnt
- Machine learning

Each of those steps add some file in the population_ directories that should have been set up at the  begininning of the procedure. Those directiories, at the beginning of the analysis, have to include:

- genomes directory:
  a directory containing the unzipped .fasta files that encode the genome sequences that will be used in order to simulate our metagenomic sample. The genomes that can be found in the example_poulation are taken from PATRIC database, like all the genomes that was used in our analysis.
- input.json:
  a .json file containing various information and metadata usefull for the simulation and the following steps of the analysis, see the example in the example_population folder
- input.tsv:
  a .tsv file containing vearious information usefull for the simulation, see the example in the example_population folder
  


## Read library generation

The read library generation is the procedure that we implemented in order to obtain simulated read libraries starting from genome sequences (.fasta files). This job is handled by Camisim and the command which is passed to it is reported in the Snakefile under the rule "camisim". In this rule are generated the reads in the anonimus mode of the software and then is used a script in order to separe the backvard and forward reads in the two pair end files. This rule is likely to require the execution of the rule "population set up". In this rule is lunched a script that generates theinput files required by Camisim, starting from the information contained in input.json and input.tsv. Following this procedure, many paramethers of Camisim execution are set with default values that we found appropriate for our goal. Anyway, if someone intend to modify those parameters should modify the script "my_scripts/generate_files.py" and follow camisim manual.

## Allignemnt

The allignment is the procedure that associates the read libraries to the references from CARD and WildCARD. This procedure is handled by the software rgi in its bwt mode. The command used to run the allignment is reported in the rule "rgi_bowtie2" and requires as input the read libraries resulting from the grneration part. Rgi requires to set up an envoirment that include the borrows weeler transform of the reference database. This procedure can be performed justo once for all, usually if is the firt time running the pipeline. If it is required, the rule "rgi_bowtie2" calls the rule "rgi_env_setup" that contains the rgi commands used in order to prepare this envoirment.

## Machine learning

This part of the pipeline is not included in the Snakefile, but is handeld by the script M_learning.py. To make this script work is important to re-organize and re-name properly the files in new directories: M_learning requires two of them, containing the same number of files. Those files should be the results of the allignment in one directory and the input.json file in the other. Those files, althought being in different directories, should be renamed and should have the same name (both are in the json format). If the set up is correct, the script performs the operations described in its documentation considering "phenotype" field of the input.json file as the ground throot about the data and the result of the allignment as the data themself. More information on this script can be found in my thesis or in the comments on the functions of the script.

The snakefiles contains also other sules, which are not already functioning, so do not realy on them.
