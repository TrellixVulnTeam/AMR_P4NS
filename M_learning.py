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
from os.path import isfile, join
from sklearn.ensemble import RandomForestClassifier
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNet
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import LeaveOneOut

average_read_length = 170

def read_bwt_2(bwt_dir):
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
	files = [f for f in os.listdir(meta_dir) if isfile(join(meta_dir, f))]
	metadata = {}
	for F in files:
		with open(join(meta_dir, F), mode = "r") as f:
			M = json.load(f)
			metadata[F.split(".")[0]] = M[list(M.keys())[0]][key]
	return metadata

def filter_sparcity_median(dataset, sparsity_filter, median_filter, plot = False, plot_zoom = 1):
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
	#FAI GRAFICO PER UN EQUA DISTRIBUZIONE DELLA PERDITA DI INFORMAZIONE
	return new

def view_database(dataset, view_a_spectrum = None):
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
	temporary = dataset.replace(to_replace = 0, value = 0.5)
	temporary = temporary.drop("phenotype", axis = 1)
	ph = dataset["phenotype"].tolist()
	G = pow(temporary.product(axis = 1, skipna = True), 1/dataset.shape[1])
	temporary = np.log(temporary.div(G, axis = 0))
	temporary["phenotype"] = ph
	return temporary

def PCA_maker(dataset, numero_componenti = None, plots = False):
	pca = PCA()	
	pca.fit(dataset.drop("phenotype", axis = 1).T)
	D = pd.DataFrame(pca.components_, columns=dataset.index.values.tolist()).T
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
	scale = StandardScaler()
	nuovo = scale.fit_transform(dataset.drop("phenotype", axis = 1))
	nuovo = pd.DataFrame(nuovo, index = list(dataset.index), columns = list(dataset.drop("phenotype", axis = 1).columns))
	nuovo.insert(nuovo.shape[1], "phenotype", dataset.loc[:,"phenotype"])
	return nuovo
	

def create_train_test(dict_di_dataframe, phenotypes, sparsity=70, median=5):
#CREAZIONE DI UN PANDAS DATAFRAME NOME (CAMPIONE X GENE RESISTENTE) CON TUTTE LE ENTRATE NECESSARIE
	ind=[]
	ll = []
	for a in dict_di_dataframe:
		ind += list(dict_di_dataframe[a].index)
	ind = list(set(ind))
	for d in dict_di_dataframe:
		new_data = {}
		for i in ind:
			if i in dict_di_dataframe[d].index:
				new_data[i] = dict_di_dataframe[d].loc[i,"FPKB"]
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
	view_database(dataset)
#FILTERING FOR PREVALENCE #toglgie quelle con occorrenza di zeri maggiore del tot%
	dataset = filter_sparcity_median(dataset, sparsity, median, plot = False, plot_zoom=15)
	view_database(dataset, 2)
#NORMALIZE #normalizza con center og rateo
	dataset = CLR_normalization(dataset)
#MAKING STANDARDIZATION
	dataset = std_scale(dataset)
#MAKING PCA #il secondo argomento standard determina quante componenti tiene, in questo caso tutte
	dataset = PCA_maker(dataset, plots=True)
	return dataset

dati = read_bwt_2("bwt")
phenotipi = read_metadata("metadata", "Phenotype")
Dati = create_train_test(dati, phenotipi)
Resistome = Dati.drop("phenotype", axis=1)
Phenotype = Dati["phenotype"]
Ground_throot_bool = [1 if d == "resistant" else 0 for d in Dati["phenotype"]]
Ground_throot = list(Dati["phenotype"])

Forest = RandomForestClassifier()
Net = ElasticNet(random_state=0)
loo = LeaveOneOut()
scores_loo_forest = cross_val_score(Forest, Resistome, Phenotype, cv=loo, n_jobs=-1)
results_loo_forest = list(cross_val_predict(Forest, Resistome, Phenotype, cv=loo, n_jobs=-1))
scores_cv5_forest = cross_val_score(Forest, Resistome, Phenotype, cv=5, n_jobs=-1)
results_cv5_forest = list(cross_val_predict(Forest, Resistome, Phenotype, cv=5, n_jobs=-1))

#scores_loo_net = cross_val_score(Net, Resistome, Sus_zero_Res_uno, cv=loo, n_jobs=-1)
#results_loo_net = cross_val_predict(Net, Resistome, Sus_zero_Res_uno, cv=loo, n_jobs=-1)
#scores_cv5_net = cross_val_score(Net, Resistome, Sus_zero_Res_uno, cv=5, n_jobs=-1)
#results_cv5_net = cross_val_predict(Net, Resistome, Sus_zero_Res_uno, cv=5, n_jobs=-1)