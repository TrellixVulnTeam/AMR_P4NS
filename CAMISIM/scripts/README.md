# CAMISIM (new) 'known distribution' modality

## Why this new modality
The modality called **known distribution** is built starting from the _de novo_ mode of [CAMISIM](https://github.com/CAMI-challenge/CAMISIM).
In the **known distribution** modality, as well as in all the other already available modalities in CAMISIM, new strains can be generated through sgEvolver.

What this new modality does differently is the community design step, which, in this case, is based on a distribution given by the user, and not randomly generated 
from a log-normal distribution.

Once the new strains are generated, this new modality will distribute the relative abundances of each "original" genome to all its simulated strains. An example of this process 
can be found in the [Wiki](https://github.com/Ettore1024/Metagenomic_known_distribution/wiki#redistribution-of-abundances).

# Installation
In this repository you can find the "_scripts_" folder of the CAMISIM tool, changed in a way that allows you to have the above-described extra modality in the _de novo_ mode.

Once you have installed CAMISIM (following the instructions here: [CAMISIM installation](https://github.com/CAMI-challenge/CAMISIM/wiki/User-manual#installation)), 
to enable this modality you may simply launch the following commands (**from inside the CAMISIM folder** created during the tool installation):

    rm -r scripts
    git clone https://github.com/Ettore1024/Metagenomic_known_distribution.git
    mv Metagenomic_known_distribution scripts

# Documentation
A complete documentation for the original CAMISIM tool can be found [here](https://github.com/CAMI-challenge/CAMISIM/wiki/User-manual).

The following sections will only be about the new **known distribution** modality and its features and options.

## What to pay attention to
For using this new modality, a new file, and a new parameter in the configuration file are required.
Also, the (already existing) _mode_ parameter should be set equal to _known_distribution_, which is the new label representing this modality.

An example of the new configuration file can be found in the [Wiki](https://github.com/Ettore1024/Metagenomic_known_distribution/wiki#example-of-configuration-file).

### The new _abundance.tsv_ file
The new file must be a _.tsv_ file with no header and two columns: the first one with the _genome_ID_ used in the other input files required by CAMISIM 
(_metadata.tsv_ and _genome_to_id.tsv_); the second one with the relative abundance of each "original" genome.

The parameter _num_real_genomes_ can of course be set to a number smaller than the number of genomes that are available in input. In this case, the _abundance.tsv_ file can be filled in two different ways:

  1. You can only list the genomes you want to use, with their relative abundance, omitting all the other genomes;

  2. You can list all the genomes you have in input; in this case, all the genomes you do not want to use during the simulation must be put at the bottom of the list 
and their relative abundance should be set to 0. 

### The new configuration file parameter 
The new parameter _path_to_abundance_file_ of the configuration file should be set equal to the absolute path of the above-mentioned _abundance.tsv_ file.
When you don't need to work in the **known distribution** modality, this file is not required, so that you can just leave the parameter blank.

## Other minor differences from the original CAMISIM tool
There are some minor differences between the original CAMISIM tool and this new "version".

First, the order of the genomes in the metadata file and in the abundance file must be the same (at least the shared set of genomes).

Second, also for the other (original) modalities, the user must pay attention to the order of the genomes in the metadata file. In fact, when the parameter _num_real_genomes_ is set to a number smaller than the number of genomes in the metadata file, the genomes considered during the simulation will only be the first _k_ ones, with _k_ equal to the value of _num_real_genomes_. 
