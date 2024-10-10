#!/home/lynden/anaconda3/envs/develop/bin/python

# Begining of Python
import pandas as pd
import random
import csv
import data_extraction
import os

# Running globally since the dataframe is required for function of script

# This is reading the csv file that contains all of the SWISS-PROT proteins which 
# is made by the data_extraction.py script
try:
    proteindata = pd.read_csv("data.csv")
except: 
    try:
    # Replace an old possilbe corrupted version if one exists.
        os.rename("data.csv", "data.csv.old")
    except:
        pass
    # Regathers the data whether that be from the uniprot_sprot.dat file or 
    # downloading that from the web again.
    print("Data was not found. Recollecting now.")
    data_extraction.get_data(outfile="data.csv")
    proteindata = pd.read_csv("data.csv")


# list contents of pd.dataframe return as strings, this uses some simiply text 
# parsing to turn it into the correct thing. 
def parse_string_to_list(mystring:str) -> list:
    return (mystring[1:-1].replace("'","").split(","))

# Creates an Row object to similfy some processing later on.
class Row():
    def __init__(self, no):
        self.content = proteindata.iloc[no]
        self.taxonomy = parse_string_to_list(self.content["taxonomy"])
        self.refs = parse_string_to_list(self.content["db-references"])
        self.no_ref = len(self.refs)
        self.kingdom = self.taxonomy[0]


# Get a selection of proteins.
def get_proteins(proteindata,number: int = 1000,GOterm_Label_No: int = 1, csv_out:str="output.csv")-> list:
    len_rows = len(proteindata.index)-1
    
    # Keeps track of the proteins that are being used in a memory efficient manner.
    protein_index = []
    
    # Counter for different kingdoms
    distribution = {"Archaea":0,"Bacteria":0,"Eukaryota":0,"Viruses":0}
    
    # Randomly selects indexes from the dataframe and then tests if there is a GO Annotation, 
    # if so then it will append it to the list else it attempts again, it also makes sure 
    # that there are no duplicates.
    for x in range(number):
        # Could be a more efficent method but it works.
        no = random.randint(0, len_rows) 
        # Each row of the dataframe
        row = Row(no)
        # conditions that must be met to be considered an eligble protein.
        checks = [distribution[row.kingdom]>=(number/4),no in protein_index,row.no_ref <= GOterm_Label_No]
        # Repeats until conditions are met otherwise skipped if they are.
        while any(checks):
            no = random.randint(0, len_rows)
            row = Row(no)  
            checks = [distribution[row.kingdom]>=(number/4),no in protein_index, row.no_ref <= GOterm_Label_No]
        
        # Add to distribution counter
        distribution[row.kingdom] += 1
        protein_index.append(no)
    # Checks the CSV file out parameter, if a non-zero length str is added as a 
    # parameter then it will write the data from the output to a csv file.
    if len(csv_out) > 0:
        with open(csv_out, "w") as f:
            csvwriter = csv.writer(f)
            # Takes already existing columns of last proteindata file and uses those.
            csvwriter.writerow(proteindata.columns)
            # Interates over the list of proteins and writes them
            for x in protein_index:
                csvwriter.writerow(proteindata.iloc[x].values)
    else: 
        # Something a dumb person does...
        print("Dumbass a zero length string for a filename isn't possible")
        raise SystemExit(1)
    # Returns the list of indexes if at some stage needed.
    return protein_index

# Used only for testing purposes to ensure that the information is correct and 
# there weren't any computation errors
def validation(path:str):
    outputdata = pd.read_csv(path,header=0)
    distribution = {"Archaea":0,"Bacteria":0,"Eukaryota":0,"Viruses":0}
    for x in outputdata["taxonomy"]:
        kingdom = parse_string_to_list(x)[0]
        distribution[kingdom] += 1
    return distribution, len(outputdata)

# For when running as a script rather than module.
def main():
    import argparse
    import time
    # parser for the command.
    parser = argparse.ArgumentParser(description="Selects a number of proteins from the SwissProt database")
    parser.add_argument("-n", "--number", type=int, help="The number of proteins to select", default=1000)
    parser.add_argument("-g", "--goterms", type=int, help="Minimum Number of GOTerms allowed", default=1)
    parser.add_argument("-o", "--output", type=str, help="Output CSV file", default="output.csv")
    args = parser.parse_args()
    # Function to run the protein collection from the dataset
    start = time.time()
    print("Processing Has Begun")
    get_proteins(proteindata,args.number, args.goterms, args.output)
    end = time.time()
    print(f"Completed in {end-start} seconds")

if __name__ == "__main__":
    main()


