#!/home/lynden/anaconda3/envs/develop/bin/python

import pandas as pd
import os
import math
import csv
import GOToolsNeo4j

# For readablity it has the following function
# Its purpose is simiply to convert a string structured as a list into a list
def parsestringintolist(mystring:str) -> list: 
    return (mystring[1:-1].replace("'","").split(","))

# File naming for fasta files
def FastaName(prefix:str,no:int):
    return f"./fasta/{prefix+str(no)}.fasta"

# Creates a object that is called ProteinData, which represents the dataset and
# thus allowing the user to perform a set of functions on the dataset.
class ProteinData:
    def __init__ (self, file_path: str):
        global df
        df = pd.read_csv(file_path,memory_map=True)
        self.df = pd.read_csv(file_path,memory_map=True)
        self.terms = len(self.df)
    
    # Writes a fasta file for the Dataset allowing for easy insertion into 
    # models, as they almost all take fasta files as a input.
    def write_fasta(self,prefix:str,number:int=0,max_per_file:int=0, incompatiable_characters:str | list[str]="", maximumlength:int=0) -> list:
        if number <= 0 or number > self.terms:
            number = self.terms

        if max_per_file <= 0 or max_per_file > number:
            max_per_file = number

        if maximumlength <= 0:
            def checklength(*args):
                return False
        else:
            def checklength(sequence_len):
                return sequence_len > maximumlength 


        no_files = math.ceil(number/max_per_file)
        filenames = [FastaName(prefix, x) for x in range(no_files)] 
        incompatiable = []
        for i,file in enumerate(filenames):
            try:
                os.remove(file)
            except:
                pass

            begining = i * max_per_file
            end = (i+1) * max_per_file
            if end > number:
                end = number
            
            with open(file, "a") as f:
                for protein in range(begining,end,1):
                    data = self[protein]
                    if any(x in data.sequence for x in incompatiable_characters) or checklength(data.sequence_len):
                        incompatiable.append(data.accession)
                    else:
                        f.write(self[protein].fasta())
        
        # Creates a file that doesn't met the model requirements and hence will be counted as zero.
        with open(f"fasta/{prefix}_incompatiable.csv", "w") as incompatiablefile:
            csvwriter = csv.writer(incompatiablefile)
            print(incompatiable)
            for x in incompatiable:
                csvwriter.writerow([x])
        
        return incompatiable

    # Creates a entry object in the dataset which can be returned by a later function.
    # This object has relvent attributes for all the proteins which can be 
    # gotton from the dataset.
    class Entry:
        def __init__(self,df,id:int):
            # Collection of important data of each Protein Entry
            self.name = df.iloc[id]["name"]
            # Name of protein as listed in SwissProt database
            self.accession = parsestringintolist(df.iloc[id]['accession'])[0] 
            # Parse the string as a list and for simiplity use the first one always
            self.organism = parsestringintolist(df.iloc[id]['organism'])
            # Get the organism of origin for the protein
            self.GOTerms = [x.strip() for x in parsestringintolist(df.iloc[id]['db-references'])]
            # Gets the leaf GO term
            #print(self.leafGOTerms)
            #tempGOTerms = [] 
            #for x in self.leafGOTerms:
            #    tempGOTerms.extend(GOToolsNeo4j.GOTerm(x).Ancestors())
            #self.GOTerms = list(set(tempGOTerms))
            # Gets all the GO Terms for the protein
            self.sequence = df.iloc[id]['sequence']
            # Get the sequence of the protein
            self.sequence_len = len(self.sequence)
            # Gets the sequence length of the protein
            self.taxonomy = parsestringintolist(df.iloc[id]['taxonomy'])
            # Gets that taxonomy of the organim of origin.
        
        # Returns the fasta format for an entry
        def fasta(self):
            return f">{self.accession} \n{self.sequence}\n"
        
        # Pretty prints the entry
        def __str__(self):
            return f"Name: {self.name}\nAccession: {self.accession}\nOrgansim of Origin: {self.organism}\nTaxonomy: {self.taxonomy}\nSequence: {self.sequence}\nSequence Length: {self.sequence_len}\nGO Annotations: {self.GOTerms}\n"
    
    # Gets the entry in the dataset by its id , which its row in the csv file/dataframe
    def get_entry_by_id(self, id:int):
        if id < 0 or id >= self.terms:
            raise ValueError("ID out of range")
        else:
            # I have replaced this function with another which in theory is better
            # and takes advantage of native python features, this function is being
            # keep for compatability purposes.
            #print("deprecation warning has been replaced by __getitem__ method")
            return self.Entry(df=self.df, id=id)
    
    # Requires further optomisation.
    def get_id_by_accession(self, accession:str):
        # I dislike this one but it works
        return self.df[self.df['accession'].str.contains(accession, na = False)].index[0]

    # Allows the code to interate over the dataset
    def __iter__(self):
        self.iter = 0
        return self
    
    def __next__(self):
        if self.iter < self.terms:
            self.iter += 1
            return self.get_entry_by_id(id=self.iter-1)
        else:
            raise StopIteration
    
    # Magic method to allow for indexing of the dataset
    def __getitem__(self, items):
        if isinstance(items, int):
            if items < 0 or items >= self.terms:
                raise ValueError("ID out of range")
            else:
                return self.Entry(df=self.df, id=items)
        # If given slice object
        elif isinstance(items, slice):
            # Dealing with None values
            if items.start is None:
                items = slice(0, items.stop, items.step)
            if items.stop is None:
                items = slice(items.start, self.terms, items.step)
            if items.step is None:
                items = slice(items.start, items.stop, 1)
            # Checks for out of range values
            if items.start < 0 or items.start >= self.terms or items.stop < 0 or items.stop >= self.terms:
                raise ValueError("Slice out of range")
            else:
                # Returns list of entries
                return [self.Entry(df=self.df, id=x) for x in range(items.start, items.stop, items.step)]
        # Get by accession
        elif isinstance(items, str):
            idx = self.get_id_by_accession(items)
            return self.Entry(df=self.df, id=idx)


def main():
    import argparse
    # parser for the command
    parser = argparse.ArgumentParser(description="The primary function in script mode is to create fasta files.")
    parser.add_argument("-c","--characters", type=str, help="The incompatible characters in a sequence. Give them in one string eg. 'uga'", default="")
    parser.add_argument("-a","--maxseqlength", type=int, help="The maximum sequence length.", default="0")
    parser.add_argument("-i","--input", type=str, help="Input CSV data file", default="output.csv")
    parser.add_argument("-n","--number", type=int, help="Number of proteins to be written to fasta file", default=0)
    parser.add_argument("-m","--max", type=int, help="Maximum number of proteins per file.", default=0)
    parser.add_argument("-o","--outfile", type=str, help="the prefix for the outfile names", default="output")
    args = parser.parse_args()
    data = ProteinData(args.input)
    data.write_fasta(args.outfile,args.number,args.max,list(args.characters),args.maxseqlength)


# For when being run as a script to create fasta files. 
if __name__ == "__main__":
    data = ProteinData("output.csv")
    print(data["A5A6H9"].GOTerms)
    #main()
