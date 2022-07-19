import pandas as pd
import sys

with open(str(sys.argv[1])) as contigs:
	line = contigs.readline()
	lista = []
	print("si parte")
	while line != "":
		if line[0] == ">":
			l = line.split("_")
			coverage = float(l[-1])
			lunghezza = int(l[-3])
			nome = ""
			for i in l[0:-4]:
				nome += i + "_"
			lista.append([nome,coverage,lunghezza])
		line = contigs.readline()
	dati = pd.DataFrame(lista, columns = ["name", "coverage", "lungth"])
	print(dati)
	dati.to_csv(str(sys.argv[2]), sep="\t")
