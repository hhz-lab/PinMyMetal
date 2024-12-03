import pandas as pd
import numpy as np
import os, sys
import getopt
sys.path.append(os.path.join(os.path.dirname(__file__), "../utils"))
import db_utils

# Get the directory where the script file resides
script_directory = os.path.dirname(os.path.abspath(__file__))

options, remainder = getopt.getopt(sys.argv[1:], 'i:o', ['input='])
for opt, arg in options:
    if opt in ('-i', '--input'):
        inputid = arg
pdbid=inputid
DBname="metal"+str(pdbid)

output_data = os.path.join(script_directory, '../', f"{pdbid}_nb_result")

# Get a database connection
conn = db_utils.create_connection(DBname)
cur = conn.cursor()

def fetch_data_and_process(sql_query, bins, n_bins, n_classes, output_file):
    # Execute the SQL query
    cur.execute(sql_query)
    data = cur.fetchall()

    # Convert the fetched data to a DataFrame
    df = pd.DataFrame(data, columns=['id', 'chem_class', 'distance'])

    # Function to create histogram features
    def create_histogram(group, bins, n_bins, n_classes):
        histogram = np.zeros((n_bins, n_classes))
        for _, row in group.iterrows():
            bin_index = np.digitize(row['distance'], bins) - 1
            if 0 <= bin_index < n_bins:
                histogram[bin_index, int(row['chem_class'])] += 1
                # histogram[bin_index, int(row['chem_class']) - 1] += 1
        return histogram.flatten()

    # Group by 'id' and compute histogram features
    histogram_features = df.groupby(['id']).apply(lambda group: create_histogram(group, bins, n_bins, n_classes))

    # Convert to DataFrame
    histogram_features_df = pd.DataFrame(histogram_features.tolist(), index=histogram_features.index)

    # Add 'id' column back to the DataFrame
    histogram_features_df.reset_index(inplace=True)

    # Add column names
    column_names = ['id']
    for i in range(1, n_bins + 1):
        for j in range(0, n_classes):
            column_names.append(f'bin{i}_class{j}')

    histogram_features_df.columns = column_names

    # Save the processed features to a CSV file
    histogram_features_df.to_csv(output_file, index=False)

# Define the range of bins excluding 0-2A
bins = np.linspace(2, 5, 4)  # Split the 2-5 A range into 3 bins
n_bins = len(bins) - 1
n_classes = 8  # Number of chemical classes

 
# pre sites for TEMSP
sql1 = "select id, chem_class, distance from chemtype_pre_chedh where sitetype='edh'"
sql2 = "select id, chem_class, distance from chemtype_pre_chedh where sitetype='ch'"

# Output file paths
output_file1 = os.path.join(output_data, 'chem_edh.csv')
output_file2 = os.path.join(output_data, 'chem_ch.csv')

# Fetch data, process, and save to CSV files
fetch_data_and_process(sql1, bins, n_bins, n_classes, output_file1)
fetch_data_and_process(sql2, bins, n_bins, n_classes, output_file2)
cur.close()
conn.close()
