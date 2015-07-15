import csv

fout = open('wifes_targets.txt')

with open('rvs.txt') as f:
	reader = csv.reader(f)
	writer = csv.writer(fout)
	for row in reader:
		if row[5]-row[6] < -15 or row[5]+row[6]>5:
			writer.writerow(row)

f.close()			
fout.close()
