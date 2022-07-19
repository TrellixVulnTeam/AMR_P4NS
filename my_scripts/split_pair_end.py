import sys

def adjust_descriptor(old):
	num = old.split("R")[1]
	n = num.split("/")
	new = "@Read_."+n[0]+" "+n[0]+"/"+n[1]
	return new

with open(str(sys.argv[1]), "r") as mixed, open(str(sys.argv[2])+"reads_1.fastq", "w") as r1, open(str(sys.argv[2])+"reads_2.fastq", "w") as r2:
	line = mixed.readline()
	while line != "":
		if list(line)[-2] == "1":
			r1.write(adjust_descriptor(line))
			line = mixed.readline()
			r1.write(line)
			line = mixed.readline()
			r1.write(line)
			line = mixed.readline()
			r1.write(line)
		elif list(line)[-2] == "2":
			r2.write(adjust_descriptor(line))
			line = mixed.readline()
			r2.write(line)
			line = mixed.readline()
			r2.write(line)
			line = mixed.readline()
			r2.write(line)
		else:
			print("qualcosa è andato storto_____" + line )
		line = mixed.readline()
