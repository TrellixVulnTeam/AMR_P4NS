#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 18 14:57:20 2022

@author: utente
"""

import os
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statistics
from os.path import isfile, join
from sklearn.ensemble import RandomForestClassifier
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import cross_validate
from sklearn.model_selection import LeaveOneOut
from sklearn.ensemble import AdaBoostClassifier
from sklearn.model_selection import train_test_split

average_read_length = 170

def read_kmer(kmer_dir):
	"""This function reads the output of kmc_count
	The rest of the analysis is not already supported for this kind of data
	"""
	files = [f for f in os.listdir(kmer_dir) if isfile(join(kmer_dir, f))]
	kmer_collection_of_samples = {}
	for F in files:
		with open(join(kmer_dir, F), mode = "r") as f:
			kmer_collection_of_samples = pd.read_csv(f, sep='\t')
	return kmer_collection_of_samples
	

def read_bwt_2(bwt_dir):
	"""Read the output of rgi
	Paramethers:
		bwt_dir: the path to the directory where the files to be read are stored
	Returns:
		a two dimensional dictionary: at first level we have a sequence of sample_name:resistome
		each resistome have the structure ARO_identity:FPKB_estimated_abundance
	"""
	files = [f for f in os.listdir(bwt_dir) if isfile(join(bwt_dir, f))]
	resistomes = {}
	for F in files:
		with open(join(bwt_dir, F), mode = "r") as f:
			dati_ = json.load(f)
			dati = [di for di in dati_ if isinstance(di, dict)]
			R={}
			for aro in dati:
				if str(aro["aro_accession"]) in R.keys():
					R[str(aro["aro_accession"])] += int(aro["reads"]["all"]) * average_read_length / int(aro["reference"]["sequence_length"])
				else:
					R[str(aro["aro_accession"])] = int(aro["reads"]["all"]) * average_read_length / int(aro["reference"]["sequence_length"])
			resistomes[F.split(".")[0]] = pd.DataFrame.from_dict(R, orient='index', columns=["FPKB"])
			print(join(bwt_dir, F) + " ha saltato " + str(len(dati_)-len(dati)) + " entrate su " + str(len(dati_)))
	return resistomes

def read_metadata(meta_dir, key):
	"""Read the metadata file
	Paramethers:
		meta_dir:
			The path to the directory where the metadata are stored.
			Notice that each metadata file have to be named exaccly equal to the data file (from rgi_bwt or kmc_count) 
		key:
			This key is refered to the structure of the json file used as metadata.
			Specify which value to consider as a label in the machine learning procedure that will follow, "phenotype" in our case
	Returns:
		A dictionary with the structure ARO_identity: phenotype
	"""
	files = [f for f in os.listdir(meta_dir) if isfile(join(meta_dir, f))]
	metadata = {}
	for F in files:
		with open(join(meta_dir, F), mode = "r") as f:
			M = json.load(f)
			metadata[F.split(".")[0]] = M[list(M.keys())[0]][key]
	return metadata

def filter_sparcity_median(dataset, sparsity_filter, median_filter, plot = False, plot_zoom = 1):
	"""Perform the filter operations on the database
	The values for the median and the sparsity filters are the two positional arguments, while the standard arguments regulates the plots of the function
	See my thesis for a bether explenation of the effects of those filters
	The original dataset is not modified.
	"""
	new = dataset.copy()
	if plot:
		y = dataset.drop("phenotype", axis=1).median(axis=0)
		y = [i+0.1 if i < 0.1 else i for i in list(y)]
		x = [(c_data == 0).astype(int).sum()*100/(new.shape[0]-1) for (c_name, c_data) in new.drop("phenotype", axis=1).iteritems()]
		colors = np.random.rand(3)
		somma = dataset.drop("phenotype", axis=1).sum(axis = 0)/50
		somma = [l*plot_zoom for l in list(somma)]
		plt.scatter(x, y, s=somma, c=colors, alpha=0.5)
		plt.yscale("log")
		plt.ylabel("Median")
		plt.xlabel("Percentage of sparsity")
		plt.title("Filtering results")
		plt.plot([0,100],[median_filter,median_filter], color = "red") #median
		plt.plot([sparsity_filter,sparsity_filter],[0,1500], color = "red") #sparsity
		plt.show()		
	for (c_name, c_data) in new.drop("phenotype", axis=1).iteritems():
		zeros_percent = (c_data == 0).astype(int).sum()*100/new.shape[0]
		if zeros_percent > sparsity_filter:
			new = new.drop(c_name, axis=1)
	median = new.drop("phenotype", axis=1).median(axis=0)
	median = median < median_filter
	for (c_name, c_data) in new.drop("phenotype", axis=1).iteritems():
		if median[c_name]:
			new = new.drop(c_name, axis=1)
	print("filter_prevalence got read of " + str(dataset.shape[1]-new.shape[1]) + " columns over " +str(dataset.shape[1]))
	aro_left_out = [k for k in dataset.drop("phenotype", axis = 1).columns if k not in new.columns]
	left_out = dataset[aro_left_out]
	print("È stato rimosso " + str(left_out.to_numpy().sum() / dataset.drop("phenotype", axis = 1).to_numpy().sum()) + " del contenuto dei resistomi")
	X = list(dataset.index)
	Yprima = list(dataset.drop("phenotype", axis=1).sum(axis = 1))
	Ydopo = list(new.drop("phenotype", axis=1).sum(axis = 1))
	X_axis = np.arange(len(X))
	plt.bar(X_axis - 0.2, Yprima, 0.4, label = 'Before')
	plt.bar(X_axis + 0.2, Ydopo, 0.4, label = 'After')
	plt.xticks(X_axis, X)
	plt.xlabel("Populations")
	plt.ylabel("Spectrum integral for the population")
	plt.legend()
	plt.show()	
	return new

def view_database(dataset, view_a_spectrum = None):
	"""Produces various plots of the dataset provided as argument and also of a specific sample
	if "view_a_spectrum" have a reasonable integer value
	"""
	V = dataset.drop("phenotype", axis = 1).sum(axis = 1)
	dic = {a : V[a] for a in V.index.tolist()}
	dic = dict(sorted(dic.items(), key=lambda item: item[1]))
	plt.bar(dic.keys(),dic.values())
	plt.xlabel("Samples")
	plt.ylabel("Resistome integral (FPKB)")
	plt.xticks([])
	plt.show()
	if isinstance(view_a_spectrum, int):
		V = dataset.drop("phenotype", axis = 1)
		dic_2 = V.iloc[view_a_spectrum,:].to_dict()
		dic_2 = dict(sorted(dic_2.items(), key=lambda item: item[1]))
		plt.bar(dic_2.keys(),dic_2.values())
		plt.xlabel("ARO indentities")
		plt.ylabel("Abundance (FPKB)")
		plt.xticks([])
		plt.show()
		
def CLR_normalization(dataset):
	"""Performs the CLR normalization on the dataset in argument
	The original dataset is not modified	
	"""
	temporary = dataset.replace(to_replace = 0, value = 0.5)
	temporary = temporary.drop("phenotype", axis = 1)
	ph = dataset["phenotype"].tolist()
	G = pow(temporary.product(axis = 1, skipna = True), 1/dataset.shape[1])
	temporary = np.log(temporary.div(G, axis = 0))
	temporary["phenotype"] = ph
	return temporary

def PCA_maker(dataset, numero_componenti = None, plots = False):
	"""Performs the PCA on the dataset
	It can be set a maximal number of components to be kept, otherwise we will end up with a number of them equal to the number of samples
	If plot is set thrue the function produces also some usefull plots
	The original dataset is not modified
	"""
	pca = PCA()	
	pca.fit(dataset.drop("phenotype", axis = 1).T)
	D = pd.DataFrame(pca.components_, columns=dataset.index.values.tolist()).T
	D.columns = ["component "+str(i) for i in range(D.shape[1])]
	if plots:
		V = pd.DataFrame(pca.explained_variance_ratio_)
		a = [str(v) for v in list(V.index)]
		b = list(V[0])
		plt.bar(a,b)
		plt.xlabel("PCA components")
		plt.ylabel("Explained variance rateo")
		plt.show()	
		plt.scatter(D.iloc[:,0], D.iloc[:,1], color = ["red" if pheno == "resistant" else "blue" for pheno in list(dataset["phenotype"])])
		plt.xlabel("First component")
		plt.ylabel("Second component")
		plt.show()
	D["phenotype"]=dataset["phenotype"].tolist()
	return D
	
def std_scale(dataset):
	"""Performs the standard scaling of the valueds in the dataframe
	The original dataset is not modified
	"""
	scale = StandardScaler()
	nuovo = scale.fit_transform(dataset.drop("phenotype", axis = 1))
	nuovo = pd.DataFrame(nuovo, index = list(dataset.index), columns = list(dataset.drop("phenotype", axis = 1).columns))
	nuovo.insert(nuovo.shape[1], "phenotype", dataset.loc[:,"phenotype"])
	return nuovo
	

def create_train_test(dict_dict, phenotypes, sparsity=70, median=5):
	"""Creates the database in a format usefull for machine learning
	Paramethers:
		dict_dict: a dictionary containing other dictionary, this argument is supposed to be the returning value from read_bwt
		phenotypes: a dictionary continining the ground throoth about sample's labels, labelled with the sample names
		The standard agruments are passed to the functino filter_median_sparsity, see the comment of it
	Returns:
		A pandas dataframe containing all the components of the data (samples) after various operations are performed: we have
		-Filtering the columns for their sparsity and median vales
		-Normalization with center log rateo
		-Standard scaling
		-PCA
	Each of the steps of this function are independent and can be commented it is required not to have them
	"""
	ind=[]
	ll = []
	for a in dict_dict:
		ind += list(dict_dict[a].index)
	ind = list(set(ind))
	for d in dict_dict:
		new_data = {}
		for i in ind:
			if i in dict_dict[d].index:
				new_data[i] = dict_dict[d].loc[i,"FPKB"]
			else:
				new_data[i] = 0
		l = []
		l.append(d)
		for i in ind:
			l.append(new_data[i])
		l.append(phenotypes[d])
		ll.append(l)
	dataset = pd.DataFrame(ll, columns=["population"]+ind+["phenotype"])
	dataset.set_index("population", inplace=True)
#	view_database(dataset)
#FILTERING FOR PREVALENCE #toglgie quelle con occorrenza di zeri maggiore del tot%
	dataset = filter_sparcity_median(dataset, sparsity, median, plot = True, plot_zoom=15)
#	view_database(dataset, 2)
#NORMALIZE #normalizza con center og rateo
	dataset = CLR_normalization(dataset)
#MAKING STANDARDIZATION
	dataset = std_scale(dataset)
#MAKING PCA #il secondo argomento standard determina quante componenti tiene, in questo caso tutte
	dataset = PCA_maker(dataset, plots=False)
	return dataset

def loo(classifier, resistome, phenotype, plot = False, linear = False):
	scores = cross_val_score(classifier, resistome, phenotype, cv=LeaveOneOut(), n_jobs=-1)
	results = list(cross_val_predict(classifier, resistome, phenotype, cv=LeaveOneOut(), n_jobs=-1))
	estimators = cross_validate(classifier, resistome, phenotype, cv=LeaveOneOut(), n_jobs=-1, return_estimator=True)
	if linear:
		feat_imps = [i.coef_[0] for i in estimators["estimator"]]
	else:
		feat_imps = [i.feature_importances_ for i in estimators["estimator"]]
	feat_imps = pd.DataFrame(feat_imps, columns=resistome.columns)
	feat_imps = {"Feature importance" : feat_imps.mean(axis=0),"Standard deviation": feat_imps.std(axis=0)}
	if plot:
		plt.hist(scores)
		plt.xlabel("Scores")
		plt.ylabel("Abundance")
		plt.show()
		X_, Y_ = [i.split(" ")[1] for i in feat_imps["Feature importance"].index], feat_imps["Feature importance"].values
		if linear:
			Y_ = [abs(iii) for iii in Y_]
		plt.bar(X_, Y_)
		plt.xlabel('PCA component')
		plt.ylabel('Feature importance')
		plt.errorbar(X_, Y_, feat_imps["Standard deviation"].values, fmt=' ',color='Black')
		plt.show()
	false_s , false_r = 0,0
	for a in range(len(phenotype)):
		if phenotype[a] == "resistant" and results[a] == "susceptible":
			false_s +=1
		elif phenotype[a] == "susceptible" and results[a] == "resistant":
			false_r +=1
	print("Leave one out report:")
	print("False resistant = " + str(false_r))
	print("False susceptible = " + str(false_s))
	print("Scores = " + str(scores))
	print("Score avrage " + str(sum(scores)/len(scores)) + ", variance " + str(np.var(scores)) + '\n')
	return {"Scores":scores, "Feature importance" : feat_imps, "Falsi suscettibili":false_s, "Falsi resistenti":false_r}

def cv5(classifier, resistome, phenotype, plot = False, linear = False):
	scores = cross_val_score(classifier, resistome, phenotype, cv=5, n_jobs=-1)
	results = list(cross_val_predict(classifier, resistome, phenotype, cv=5, n_jobs=-1))
	estimators = cross_validate(classifier, resistome, phenotype, cv=5, n_jobs=-1, return_estimator=True)
	if linear:
		feat_imps = [i.coef_[0] for i in estimators["estimator"]]
	else:
		feat_imps = [i.feature_importances_ for i in estimators["estimator"]]
	feat_imps = pd.DataFrame(feat_imps, columns=resistome.columns)
	feat_imps = {"Feature importance" : feat_imps.mean(axis=0),"Standard deviation": feat_imps.std(axis=0)}
	if plot:
		plt.hist(scores)
		plt.xlabel("Scores")
		plt.ylabel("Abundance")
		plt.show()
		X_, Y_ = [i.split(" ")[1] for i in feat_imps["Feature importance"].index], feat_imps["Feature importance"].values
		if linear:
			Y_ = [abs(iii) for iii in Y_]
		plt.bar(X_, Y_)
		plt.xlabel('PCA component')
		plt.ylabel('Feature importance')
		plt.errorbar(X_, Y_, feat_imps["Standard deviation"].values, fmt=' ',color='Black')
		plt.show()
	false_s , false_r = 0,0
	for a in range(len(phenotype)):
		if phenotype[a] == "resistant" and results[a] == "susceptible":
			false_s +=1
		elif phenotype[a] == "susceptible" and results[a] == "resistant":
			false_r +=1
	print("5 fold report:")
	print("False resistant = " + str(false_r))
	print("False susceptible = " + str(false_s))
	print("Scores = " + str(scores))
	print("Score avrage " + str(sum(scores)/len(scores)) + ", variance " + str(np.var(scores)) + '\n')
	return {"Scores":scores, "Feature importance" : feat_imps, "Falsi suscettibili":false_s, "Falsi resistenti":false_r}
	

dati = read_bwt_2("bwt_plot")
#dati = read_bwt_2("bwt_4d")
phenotipi = read_metadata("metadata_plot", "Phenotype")
#phenotipi = read_metadata("metadata_4d", "Phenotype")
Dati = create_train_test(dati, phenotipi)
Resistome = Dati.drop("phenotype", axis=1)
Phenotype = Dati["phenotype"]
Phenotype_bool = pd.Series({k : (1 if Phenotype[k] == "resistant" else 0) for k in Phenotype.keys()})

print("\nForest classifier:")
Forest = RandomForestClassifier(random_state=0)
Forest_loo_scores = loo(Forest, Resistome, Phenotype, plot=True)
Forest_cv5_scores = cv5(Forest, Resistome, Phenotype, plot=True)

print("\nAdaBoost classifier:")
Ada = AdaBoostClassifier(n_estimators=100, random_state=0)
Ada_loo_scores = loo(Ada, Resistome, Phenotype, plot=True)
Ada_cv5_scores = cv5(Ada, Resistome, Phenotype, plot=True)

print("\nLog regression classifier:")
Log = LogisticRegression(random_state=20, penalty="elasticnet", solver = "saga", l1_ratio=0.5)
Net_loo_scores = loo(Log, Resistome, Phenotype, plot=True, linear=True)
Net_cv5_scores = cv5(Log, Resistome, Phenotype, plot=True, linear=True)