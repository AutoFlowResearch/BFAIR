#CONVERT A CSV TO A DICTIONARY IN PYTHON
import csv
import sys
import pprint

# Function to convert a csv file to a list of dictionaries.
def csv_dict_list(variables_file):
     
    # Open variable-based csv, iterate over the rows and map values to a list of dictionaries containing key/value pairs
 
    reader = csv.DictReader(open(variables_file, 'rb'))
    dict_list = []
    for line in reader:
        dict_list.append(line)
    return dict_list
 
# Calls the csv_dict_list function, passing the named csv
 
device_values = csv_dict_list(sys.argv[1])
 
# Prints the results
 
pprint.pprint(device_values)

# Convert list to dict to make sure 

def Convert(device_values): 
	device_values_dct = {device_values[i]: device_values[i + 1] for i in range(0, len(device_values), 2)} 
	return device_values_dct 
		
#INSERT PYTHON DICT AFTER CONVERTION FROM CSV USING PSYCOPG2

# import psycopg2 
from psycopg2.extensions import AsIs

# connect the database
conn=psycopg2.connect("dbname='db_name' user='db_user' password='db_pass'")

# device_values_dct is dict (already converted)

columns = device_values_dct.keys()
values = [device_values_dct[column] for column in columns]

# insert statement to sql 
insert_statement = 'insert into device_values_dct_table (%s) values %s'

    # cursor.execute(insert_statement, (AsIs(','.join(columns)), tuple(values)))
print cursor.mogrify(insert_statement, (AsIs(','.join(columns)), tuple(values)))

conn.commit()
