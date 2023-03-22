import csv

file_data = "file_data.txt"

with open(file_data, 'r') as data_file:
    reader = csv.reader(data_file, delimiter=' ')
    for row in reader:
        import pdb
        pdb.set_trace()
    
