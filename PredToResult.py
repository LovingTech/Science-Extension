import GOToolsNeo4j
import GOSemanticSimilarity
import proteintools
import multiprocessing as mp
import os
import time
import csv

# initalises the particular csv dataset that will be used.
data = proteintools.ProteinData("output.csv")

# used for turing the file into a useable type
def SPROFchecker(pred_raw:list):
    returnlist = []
    for x in pred_raw:
        term = x.split("|")[0][:-1].strip()
        if "GO" in term:
            returnlist.append(term)
    return returnlist

# Convert the SPROF_TOP result into a useable result.
def SPROF_TOP(file:str):
    with open(file,"r") as f:
        proteins = f.read().split("\n\n")
    for i,x in enumerate(proteins):
        entry = x.split("\n")
        GoTerms = SPROFchecker(entry[2:])
        if None in GoTerms:
            print("error here")
        proteins[i] = (entry[0][:-1],GoTerms)
    return proteins

def NETGO3(filename:str):
    data = {}
    proteins = []
    os.system(f"""sed "s/=====//" -i {filename}""")
    with open(filename,"r") as file:
        protdata = csv.reader(file, delimiter="\t")
        for x in protdata:
            if len(x) == 0:
                pass
            else:
                if x[0] in data:
                        data[x[0]].append(x[1])
                else:
                    data[x[0]] = [x[1]] 
    for y in data:
        proteins.append((y,data[y]))
    return proteins

def DeepGoPlus(file:str):
    proteins = [] 
    with open(file,"r") as f:
        proteinlines = csv.reader(f, delimiter="\t")
        for x in proteinlines:
            proteins.append((x[0],[go[:-6] for go in x[1:]]))
    return proteins

# for dealing with single proteins
def CompareProteins(entry:tuple):
    accession,terms = entry
    label = StrToTermsIter(data[accession].GOTerms)
    preds = StrToTermsIter(terms)
    similarity = GOSemanticSimilarity.ProteinFunctionSimilarity(label,preds).Average()#.neo4jspeedmode()
    print(accession,"\t\t",similarity)
    return (accession, similarity,terms)

# converts a list Go terms strings into go term objects
def StrToTermsIter(result) -> list[GOToolsNeo4j.GOTerm]:
    return GOToolsNeo4j.list_to_terms(result)

# By using multi core map function, the process can be sped up. 
def compareall(items,outfile:str, threads:int=10):
    #dealing with the data
    print("Beginning Comparision")
    print("Accession \t Similarity")
    #input is of type list[tuple[str,list]] called items
    with mp.Pool(threads) as p:
        similarity = p.map(CompareProteins,items,30)
    
    try:
        os.remove(outfile)
    except:
        pass

    with open(outfile, "w") as out:
        w = csv.writer(out)
        w.writerows(similarity)
    return similarity

# TESTING
if __name__ == "__main__":
    start = time.time()
    DeepGOpl = DeepGoPlus("./modelresults/results_deepgo.tsv")
    compareall(DeepGOpl,"./results_updated_v2/results_deepgopl.csv",threads=32)
    sprof0 = SPROF_TOP("./modelresults/SPROF/SPROF0_top_preds.txt")
    compareall(sprof0,"./results_updated_v2/results_sprof0.csv",threads=32)
    sprof1 = SPROF_TOP("./modelresults/SPROF/SPROF1_top_preds.txt")
    compareall(sprof1,"./results_updated_v2/results_sprof1.csv",threads=32)
    for x in range(10):
        netgo = NETGO3(f"modelresults/netgo3/result_netgo{x}.txt")
        compareall(netgo,f"./results_updated_v2/results_netgo{x}.csv",threads=32)
    end = time.time()
    print("Time taken to compare predictions:", (end-start)/60)
