import csv
import json
import sys
import os
from configparser import ConfigParser

########################## Get inputs' directory ############################

try:
	working_input_dir = sys.argv[1]
except IndexError:
	print('''\nThe absolute path of the directory containing input.tsv and input.json was not passed.\n
Press [y] if the above-mentioned directory is the one you are working in:''')
	if input() == 'y':
		working_input_dir = os.getcwd()
	else:
		print('No directory given, the run will stop!')
		sys.exit()

########################## Read .tsv input file #############################

input_tsv_name = '/input.tsv'
input_tsv_path = working_input_dir + input_tsv_name

with open(input_tsv_path, 'r') as input_tsv:
	read_input_tsv = csv.reader(input_tsv, delimiter = '\t')
	next(read_input_tsv, None)

	genome_ID_list = []
	species_list = []
	abundance_list = []
	patric_ID_list = []
	antibiotic_list = []
	phenotype_list = []
	genome_len_list = []
	novelty_cat_list = []
	genome_filename_list = []

	for rows in read_input_tsv:
		genome_ID_list.append(rows[0])
		species_list.append(rows[1])
		abundance_list.append(rows[2])
		patric_ID_list.append(rows[3])
		antibiotic_list.append(rows[4])
		phenotype_list.append(rows[5])
		genome_len_list.append(rows[6])
		novelty_cat_list.append(rows[7])
		genome_filename_list.append(rows[8])

########################## Read .json input file ############################

input_json_name = '/input.json'
input_json_path = working_input_dir + input_json_name

with open(input_json_path, 'r') as input_json:
	json_to_dict = json.load(input_json)

simul_dir = json_to_dict['simulation_directory']

######################### Write abundance file ##############################

abundance_name = '/abundance.tsv'
abundance_path = simul_dir + abundance_name

with open(abundance_path, 'w') as abundance_file:
	wrt_abundance = csv.writer(abundance_file, delimiter='\t')
	for index in range(len(genome_ID_list)):
		wrt_abundance.writerow([genome_ID_list[index], abundance_list[index]])

########################## Write metadata file ##############################

metadata_name = '/metadata.tsv'
metadata_path = simul_dir + metadata_name

with open(metadata_path, 'w') as metadata_file:
	wrt_metadata = csv.writer(metadata_file, delimiter='\t')
	wrt_metadata.writerow(['genome_ID', 'OTU', 'NCBI_ID', 'novelty_category'] )
	for index in range(len(genome_ID_list)):
		wrt_metadata.writerow([genome_ID_list[index], str(index+1), str(int(float(patric_ID_list[index]))), novelty_cat_list[index]])

######################### Write genome_to_id file ###########################

genome_to_id_name = '/genome_to_id.tsv'
genome_to_id_path = simul_dir + genome_to_id_name
with open(genome_to_id_path, 'w') as genome_to_id_file:
	wrt_genome_to_id = csv.writer(genome_to_id_file, delimiter='\t')
	for index in range(len(genome_ID_list)):
		wrt_genome_to_id.writerow([genome_ID_list[index], str(simul_dir+'/genomes/'+genome_filename_list[index])])

######################## Write configuration file ###########################

def config_definition (configuration_file):

	configuration_file['Main'] = {
                 'seed': '',
                 'phase': '0',
                 'max_processor': '',
                 'dataset_id': 'RL',
                 'output_directory' : '',
                 'temp_directory': '/tmp',
                 'gsa': 'False',
                 'pooled_gsa': 'False',
                 'anonymous': 'True',
                 'compress': '1'

        }

	configuration_file['ReadSimulator'] = {
                 'readsim': 'CAMISIM/tools/art_illumina-2.3.6/art_illumina',
                 'error_profiles': 'CAMISIM/tools/art_illumina-2.3.6/profiles',
                 'samtools': 'CAMISIM/tools/samtools-1.3/samtools',
                 'profile': 'mbarc',
                 'size': '',
                 'type': 'art',
                 'fragments_size_mean': '',
                 'fragment_size_standard_deviation': ''
        }

	configuration_file['CommunityDesign'] = {
                 'ncbi_taxdump': 'CAMISIM/tools/ncbi-taxonomy_20170222.tar.gz',
                 'strain_simulation_template': 'CAMISIM/scripts/StrainSimulationWrapper/sgEvolver/simulation_dir',
                 'number_of_samples': '',
	}

	configuration_file['community0'] = {
                 'metadata': '',
                 'id_to_genome_file': '',
                 'id_to_gff_file': '',
                 'path_to_abundance_file': '',
                 'genomes_total': '',
                 'num_real_genomes': '',
                 'max_strains_per_otu': '1',
                 'ratio': '1',
                 'mode': 'known_distribution',
                 'equally_distributed_strains': '',
                 'input_genomes_to_zero': 'True',
                 'log_mu': '1',
                 'log_sigma': '2',
                 'gauss_mu': '1',
                 'gauss_sigma': '1',
                 'view': 'False'
	}
	return configuration_file

def config_settings(configuration_file, input_dict):
	config.set('Main', 'max_processor', input_dict['max_processor'])
	config.set('Main', 'seed', input_dict['seed'])
	config.set('Main', 'output_directory', input_dict['output_directory'])
	config.set('ReadSimulator', 'size', input_dict['size'])
	config.set('ReadSimulator', 'fragments_size_mean', input_dict['fragm_size_mean'])
	config.set('ReadSimulator', 'fragment_size_standard_deviation', input_dict['fragm_size_std_dev'])
	config.set('CommunityDesign', 'number_of_samples', input_dict['number_of_samples'])
	config.set('community0', 'metadata', metadata_path)
	config.set('community0', 'id_to_genome_file', genome_to_id_path)
	config.set('community0', 'path_to_abundance_file', abundance_path)
	config.set('community0', 'genomes_total', input_dict['genomes_total'])
	config.set('community0', 'num_real_genomes', str(len(genome_ID_list)))
	config.set('community0', 'equally_distributed_strains', input_dict['equally_distributed_strains'])

config_name = '/config.ini'
config = ConfigParser()

config_definition(config)
config_settings(config, json_to_dict)

config_path = simul_dir + config_name
with open(config_path, 'w') as config_file:
	config.write(config_file)

###################### Write genomes' informations #########################

tsv_info_to_dict = {}
for index in range(len(genome_ID_list)):
	tsv_info_to_dict[genome_ID_list[index]] = {}

with open(input_tsv_path, 'r') as input_tsv:
	input_tsv_reader = csv.DictReader(input_tsv, delimiter = '\t')
	i=0
	for rows in input_tsv_reader:
		tsv_info_to_dict[genome_ID_list[i]].update(dict(rows))
		i += 1

keys_list = ['Genome_ID', 'Novelty_category']
for index in range(len(genome_ID_list)):
	for key in keys_list:
		try:
			del tsv_info_to_dict[genome_ID_list[index]][key]
		except KeyError:
			pass

out_genomes_data_name = '/genomes_metadata.json'
out_genomes_data_path = simul_dir + out_genomes_data_name

with open(out_genomes_data_path, 'w') as out_genomes_data_file:
	json.dump(tsv_info_to_dict, out_genomes_data_file, indent = 8)
