import Bio.SwissProt as sp
import csv
import os
import termcolor as tc

def get_data(filepath:str="uniprot_sprot.dat", outfile:str="data.csv"):
    if os.path.exists(filepath):
        # Checks if a file exists and it can be parsed corretly
        file = sp.parse(filepath)
    else:
        # If the file cannot be found or parse correctly it will redownload the latest version and parse that.
        print(tc.colored("Error: File not found Downloading Now","red"))
        os.system("wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz")
        os.system("gzip -d uniprot_sprot.dat.gz")
        if os.path.exists(filepath):
            print(tc.colored("Download Complete","green"))
            file = sp.parse(filepath)
    
    # Outfile is written with name etc
    with open(outfile, 'w', newline='') as csvfile:
        datawriter = csv.writer(csvfile)
        # Writes table headers
        datawriter.writerow(["name", "accession","organism", "taxonomy","sequence", "db-references"])
        for record in file:
            go = [e[1] for e in record.cross_references if "GO" in e]
            datawriter.writerow([record.entry_name,record.accessions,record.organism,record.organism_classification,record.sequence, go])

# Testing
if __name__ == "__main__":
    get_data(outfile="testing.csv")
