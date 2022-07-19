__author__ = 'hofmann'
__original_author__ = 'eik.dahms@uni-duesseldorf.de'
__version__ = '0.0.3'

import random
import csv
import numpy as np
from scripts.Validator.validator import Validator
from datetime import datetime
#from decimal import Decimal

class PopulationDistribution(Validator):
	"""
	applylognorm
	inital draws for each entry/strain from a log-norm distribution

	"""
	_label = "PopulationDistribution"

	_modi = [
		'known_distribution',
		'replicates',
		'timeseries_normal',
		'timeseries_lognormal',
		'differential']

	def __init__(self, logfile=None, verbose=True, debug=False, seed=None):
		"""
			Initialize instance with seed

			@attention:

			@param logfile: file handler or file path to a log file
			@type logfile: str | file | io.FileIO | StringIO.StringIO
			@param verbose: Not verbose means that only warnings and errors will be past to stream
			@type verbose: bool
			@param debug: If True logger will output DEBUG messages
			@type debug: bool
			@param seed: The seed used for initiation of the 'random' module
			@type seed: long | int | float | str | unicode

			@return: None
			@rtype: None
		"""
		assert isinstance(verbose, bool)
		assert isinstance(debug, bool)
		super(PopulationDistribution, self).__init__(logfile, verbose, debug)

		if seed is not None:
			random.seed(seed)

	@staticmethod
	def get_valid_modes():
		return PopulationDistribution._modi

	@staticmethod
	def _get_initial_list(size_of_population, number_of_samples):
		"""
			Get initial list with zero initialized

			@attention: Each list in the list contains the distribution value for each sample

			@param size_of_population: Amount of genomes or individuals
			@type size_of_population: int | long
			@param number_of_samples: Number of samples
			@type number_of_samples: int | long

			@return: A list of lists.
			@rtype: list[list[float]]
		"""
		assert isinstance(size_of_population, int)
		assert isinstance(number_of_samples, int)
		return [[0.0] * number_of_samples for _ in range(size_of_population)]

	@staticmethod
	def Broken_stick_model(genome_name, genome_abundance_dict, genome_abundance, strain_num, param_b):
		if strain_num == 0:
			return genome_abundance_dict
		else:
			v = np.random.beta(1, param_b, size=strain_num)
			v[1:] *= np.cumprod(1 - v[:-1])
			list_for_genome = list(np.insert(v, 0, 1-np.cumsum(v)[-1]))
			i = 0
			for keys in genome_abundance_dict.keys():
				if str(genome_name) in str(keys):
					genome_abundance_dict[keys] = (list_for_genome[i] * genome_abundance)
					i += 1
			return genome_abundance_dict

	@staticmethod
	def Broken_stick_model_2(genome_name, genome_abundance_dict, genome_abundance, strain_num, param_b):
		if strain_num == 0:
			return genome_abundance_dict
		else:
			v = np.random.beta(1, param_b, size=strain_num-1)
			v[1:] *= np.cumprod(1 - v[:-1])
			v = np.insert(v, 0, 0)
			v = np.insert(v, len(v), 1-np.cumsum(v)[-1])
			list_for_genome = list(v)
			i = 0
			for keys in genome_abundance_dict.keys():
				if str(genome_name) in str(keys):
					genome_abundance_dict[keys] = (list_for_genome[i] * genome_abundance)
					i += 1
			return genome_abundance_dict

	def _consider_abundance_input(self, list_population, number_of_samples, abundance_file_path, list_of_genome_id, bool_input_genomes_to_zero):
		"""
                        @author: Ettore Rocchi
                        This function takes the initial distribution given by the user through the input tsv file
                        and distributes the original genomes' abundances to all their strains

                        @attention: Abundances are distributed differently in different samples

			@param list_population: Main list for all distributions, any element of the "external" list
                                                is a list of abundances of a (original or simulated) genome, in the
                                                different samples
                        @type : list[list[float]]
                        @param number_of_samples: Number of samples
                        @type : int | long
                        @param abundance_file_path: File path to the input distribution tsv file
                        @type : str
                        @param list_of_genome_id: List of genome ids of both original and simulated strains
                        @type : list
                """
		assert isinstance(list_population, list)
		assert isinstance(number_of_samples, int)
		assert isinstance(abundance_file_path, str)
		now = datetime.now()
		genome_to_strain = {}
		genome_to_abundance_original = {}
		with open(abundance_file_path, 'r') as abundance_file:
			for rows in csv.reader(abundance_file, delimiter='\t'):
				if rows[0] in list_of_genome_id:
					genome_id = rows[0]
					genome_to_abundance_original[genome_id] = float(rows[1])
					genome_to_strain[genome_id] = 0
				else:
					if not rows[1] in ('0.0', '0'):
						print('{} WARNING: {} has abundance in the abundance file different from 0.0 (= {})'.format(now.strftime("%Y-%m-%d %H:%M:%S"), rows[0], rows[1]))
						print('The abundance of the original genomes and their strains (total of {} genomes) will be normalized to 1'.format(str(len(list_of_genome_id))))
		total_genome_abundance = genome_to_abundance_original.copy()
		for genome_ID in list_of_genome_id:
			if not genome_ID in genome_to_abundance_original.keys():
				total_genome_abundance[genome_ID] = 0.0
				for genome in genome_to_abundance_original.keys():
					if str(genome) in str(genome_ID):
						genome_to_strain[genome] += 1
		if bool_input_genomes_to_zero == False:
			for sample in range(number_of_samples):
				dict = {}
				for genome_id in genome_to_abundance_original.keys():
					dict = self.Broken_stick_model(genome_id, total_genome_abundance, genome_to_abundance_original[genome_id], genome_to_strain[genome_id], 3)
				for genome in dict.keys():
					list_population[list_of_genome_id.index(genome)][sample] = dict[genome]
		else:
			for sample in range(number_of_samples):
				dict = {}
				for genome_id in genome_to_abundance_original.keys():
					dict = self.Broken_stick_model_2(genome_id, total_genome_abundance, genome_to_abundance_original[genome_id], genome_to_strain[genome_id], 3)
				for genome in dict.keys():
					list_population[list_of_genome_id.index(genome)][sample] = dict[genome]

	@staticmethod
	def _add_initial_log_distribution(list_population, mu, sigma):
		"""
			Adding a initial distribution

			@attention: Values for first sample

			@param list_population: Main list for all distributions
			@type : list[list[float]]
			@param mu: Mean
			@type mu: float
			@param sigma: standard deviation
			@type sigma: float

			@return: Nothing
			@rtype: None
		"""
		assert isinstance(list_population, list)
		assert isinstance(mu, (float, int))
		assert isinstance(sigma, (float, int))
		for index in range(len(list_population)):
			list_population[index][0] = random.lognormvariate(mu, sigma)

	def _add_replicates(self, list_population, mu, sigma):
		"""
			Adding gaussian noise to the first drawn abundances

			@attention:

			@param list_population: Main list for all distributions
			@type : list[list[float]]
			@param mu: Mean
			@type mu: float
			@param sigma: standard deviation
			@type sigma: float

			@return: Nothing
			@rtype: None
		"""
		assert isinstance(list_population, list)
		assert isinstance(mu, (float, int))
		assert isinstance(sigma, (float, int))
		for index_p in range(len(list_population)):
			initial_log_distribution = list_population[index_p][0]
			for index_i in range(len(list_population[index_p])-1):
				list_population[index_p][index_i+1] = self.lt_zero(initial_log_distribution + random.gauss(mu, sigma))

	def _add_timeseries_gauss(self, list_population, mu, sigma):
		"""
			Adding gaussian noise sequentially to the previous sample

			@attention:

			@param list_population: Main list for all distributions
			@type : list[list[float]]
			@param mu: Mean
			@type mu: float
			@param sigma: standard deviation
			@type sigma: float

			@return: Nothing
			@rtype: None
		"""
		assert isinstance(list_population, list)
		assert isinstance(mu, (float, int))
		assert isinstance(sigma, (float, int))
		for index_p in range(len(list_population)):
			for index_i in range(len(list_population[index_p])-1):
				if list_population[index_p][index_i] > 0:
					list_population[index_p][index_i+1] = self.lt_zero(list_population[index_p][index_i] + random.gauss(mu, sigma))
				else:
					# extinction
					list_population[index_p][index_i+1] = 0.0

	@staticmethod
	def _add_timeseries_lognorm(list_population, mu, sigma):
		"""
			each abundance profile is produced by
			- draw new value from lognorm distribution
			- add old and new value and divide by 2

			@attention:

			@param list_population: Main list for all distributions
			@type : list[list[float]]
			@param mu: Mean
			@type mu: float
			@param sigma: standard deviation
			@type sigma: float

			@return: Nothing
			@rtype: None
		"""
		assert isinstance(list_population, list)
		assert isinstance(mu, (float, int))
		assert isinstance(sigma, (float, int))
		for index_p in range(len(list_population)):
			for index_i in range(len(list_population[index_p])-1):
				list_population[index_p][index_i+1] = (list_population[index_p][index_i] + random.lognormvariate(mu, sigma))/2

	@staticmethod
	def _add_differential(list_population, mu, sigma):
		"""
			Abundance is drawn independently from previous lognorm distributions

			@attention:

			@param list_population: Main list for all distributions
			@type : list[list[float]]
			@param mu: Mean
			@type mu: float
			@param sigma: standard deviation
			@type sigma: float

			@return: Nothing
			@rtype: None
		"""
		assert isinstance(list_population, list)
		assert isinstance(mu, (float, int))
		assert isinstance(sigma, (float, int))
		for index_p in range(len(list_population)):
			for index_i in range(len(list_population[index_p])-1):
				list_population[index_p][index_i+1] = random.lognormvariate(mu, sigma)

	def display_figures(self, list_population):
		"""
			Display distribution

			@attention: limited to first 20

			@param list_population: Main list for all distributions
			@type list_population: list[list[float]]

			@return: Nothing
			@rtype: None
		"""
		assert isinstance(list_population, list)
		import matplotlib.pyplot as plt
		try:
			# display histogram
			plt.figure(1)
			plt.subplot(1, 2, 1)
			for i in range(len(list_population[0])):
				hist = {}
				print_hist = [[], []]
				for genome_id_dist in list_population:
					z = int(round(genome_id_dist[i]))
					if z in hist:
						hist[z] += 1
					else:
						hist[z] = 1

				for k, v in hist.items():
					print_hist[0].append(k)
					print_hist[1].append(v)
				# print(hist)
				# print(print_hist)
				plt.plot(print_hist[0], print_hist[1])
		except Exception:
			self._logger.error("Can not display distributions!")
			raise Exception("Can not display distributions!")

		plt.subplot(1, 2, 2)
		print_num = 20 if len(list_population) > 20 else len(list_population)
		for i in range(print_num):
			v = []
			for y in list_population[i]:
				v.append(y)
			plt.plot(v)
		plt.show()

	def get_lists_of_distributions(
		self, size_of_population, number_of_samples, abundance_file_path, bool_input_genomes_to_zero, list_of_genome_id, modus, log_mu, log_sigma, gauss_mu=None, gauss_sigma=None,
		view_distribution=False):
		"""
			Get list of distributions of all samples

			@attention:

			@param size_of_population: Amount of genomes or individuals
			@type size_of_population: int | long
			@param number_of_samples: Number of samples
			@type number_of_samples: int | long
			@param abundance_file_path: File path to the input distribution file
                        @type : str
                        @param list_of_genome_id: List of genome id
                        @type : list[str | unicode]
			@param modus: 'known_distribution', 'replicates', 'timeseries_normal','timeseries_lognormal', 'differential'
			@type modus: str
			@param log_mu: Mean for log
			@type log_mu: float
			@param log_sigma: standard deviation for log
			@type log_sigma: float
			@param gauss_mu: Mean for gauss
			@type gauss_mu: float
			@param gauss_sigma: standard deviation for gauss
			@type gauss_sigma: float

			@return: Main list for all distributions
			@rtype: list[list[float]]
		"""
		assert isinstance(size_of_population, int)
		assert isinstance(number_of_samples, int)
		assert isinstance(modus, str)
		assert isinstance(log_mu, (float, int))
		assert isinstance(log_sigma, (float, int))
		assert isinstance(gauss_mu, (float, int))
		assert isinstance(gauss_sigma, (float, int))
		list_of_genome_id_new = list(list_of_genome_id) #changes a dict.keys() into a list
		assert isinstance(list_of_genome_id_new, list)
		if gauss_mu is None:
			gauss_mu = 0
		if gauss_sigma is None:
			# TODO: gauss sigma needs proper dependence of log sigma
			gauss_sigma = 3 * log_sigma

		list_population = self._get_initial_list(size_of_population, number_of_samples)

		if modus != 'known_distribution':
			self._add_initial_log_distribution(list_population, log_mu, log_sigma)

		while True:

			if modus == 'known_distribution':
				self._consider_abundance_input(list_population, number_of_samples, abundance_file_path, list_of_genome_id_new, bool_input_genomes_to_zero)
			elif modus == 'replicates':
				self._add_replicates(list_population, gauss_mu, gauss_sigma)
			elif modus == 'timeseries_normal':
				self._add_timeseries_gauss(list_population, gauss_mu, gauss_sigma)
			elif modus == 'timeseries_lognormal':
				self._add_timeseries_lognorm(list_population, log_mu, log_sigma)
			elif modus == 'differential':
				self._add_differential(list_population, log_mu, log_sigma)

			if not view_distribution:
				break

			self.display_figures(list_population)
			if self.get_confirmation(message="Use distribution? [y/n]"):
				break
		self.random_distribution_to_relative_abundance(list_population)
		return list_population

	@staticmethod
	def random_distribution_to_relative_abundance(list_population, precision=10):
		"""
			Replace random distributions with relative abundances

			@attention: limited to first 20

			@param list_population: Main list for all distributions
			@type list_population: list[list[float]]
			@param precision: Precision, numbers after decimal point
			@type precision: int
		"""
		number_of_samples = len(list_population[0])
		for index_i in range(number_of_samples):
			total_abundance = 0.0
			for index_p in range(len(list_population)):
				total_abundance += list_population[index_p][index_i]
			for index_p in range(len(list_population)):
				list_population[index_p][index_i] = round(list_population[index_p][index_i] / float(total_abundance), precision)

	@staticmethod
	def lt_zero(value):
		"""
			Prevent values <= 0

			@attention:

			@param value:
			@type value: float | int | long

			@return: value if > 0, else 0.001
			@rtype: float | int | long
		"""
		if value <= 0:
			# > 0 to prevent extinction
			return 0.001
		else:
			return value

	def get_confirmation(self, message):
		"""
			Confirm something by requesting user input

			@attention:

			@param message: A yes no question.
			@type message: str | unicode

			@return: Answer to question
			@rtype: bool
		"""
		user_input = raw_input("{}\n>".format(message)).lower()
		while True:
			if self.is_boolean_state(user_input):
				return self.get_boolean_state(user_input)
			user_input = raw_input("Please type 'n' for no, or 'y' for yes:\n>").lower()
