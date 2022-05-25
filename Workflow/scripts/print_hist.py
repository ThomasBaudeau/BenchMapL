import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-deep')
import sys

def main():

	if len(sys.argv) == 2:
		fof_name = sys.argv[1]
	else:
		print("give me a file of files with .hist of ntcard in lines, with format species_kValue.hist")
		return 0

	fof = open(fof_name, 'r')

	files = []
	for filename in fof.readlines():
		files.append(filename.rstrip())
	fof.close()
	scale = [0]
	y_list = []
	y_list_names = []
	for hist_file in files:
		hf = open(hist_file, 'r')
		y_list.append([])
		y_list_names.append(hist_file.split('.')[0])
		consecutive_zeros = 0
		for val,lines in enumerate(hf.readlines()):
			if val > 1:
				number = int(lines.split('\t')[0])
				kmer = int(lines.rstrip().split('\t')[1])
				if kmer == 0:
					consecutive_zeros += 1
				if consecutive_zeros < 5:
					y_list[-1].extend(kmer*[val-2])
					
					if scale[-1] < number:
						scale.append(number)



	plt.hist(y_list, scale, label=y_list_names)
	plt.legend(loc='upper right')
	# ~ plt.show()
	plt.savefig(y_list_names[0].split("_")[0] + ".pdf", dpi=150)


if __name__ == "__main__":
    main()
