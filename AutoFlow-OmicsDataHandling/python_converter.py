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
# read and convert to dict using pandas

	#READ TSV FROM ISATAB
import pandas as pd
data = pd.read_csv (variables_file,sep='\t')
df = pd.DataFrame(data)
#print(df) check
#print(type(df)) check

#CONVERT PANDAS DATAFRAME INTO A DICT

list_of_df = df.to_dict()
print(list_of_df)
print(type(list_of_df))
		
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

# read data from SQL and create pandas dataframe
from pandas import DataFrame
import psycopg2
# connect to database
conn=psycopg2.connect("dbname='db_name' user='db_user' password='db_pass'")
cur = conn.cursor()
# get all data via select
cur.execute("SELECT * FROM tablename ) 
the_data = cur.fetchall()
colnames = [desc[0] for desc in cur.description]
# save data in pandas.dataframe
the_data_dataframe = DataFrame(the_data)
the_data_dataframe.columns = colnames
# close psycopg2
cur.close()
conn.close()
	    
	    
# in case to create csv file from * database
# connect to DB #
from pandas import DataFrame
import psycopg2
conn = psycopg2.connect("dbname='db_name' user='db_user' password='db_pass'")
cur = conn.cursor()
# copy DB to csv
cur.execute("COPY (SELECT * FROM database) TO '/tmp/filename.csv' (format CSV)")

	    # create template Isatab.py
	    
from isatools.model import *
def create_descriptor():
    """Returns a simple but complete ISA-Tab 1.0 descriptor for illustration."""

    # Create an empty Investigation object and set some values to the instance variables.
# and then change value from csv or dataframe?
