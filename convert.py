import csv
import proteintools
import pandas as pd


proteindata = proteintools.ProteinData("output.csv")

results = "results_updated/results_deepgopl.csv"
complete = "result_deepgo_complete_updated.csv"
df = pd.read_csv(results,header=None)

with open(complete,"w") as f:
    csvwriter = csv.writer(f)
    csvwriter.writerow(["accession","similarity","organsim","taxonomy","predictions","actual","sequence"])
    for i,x in enumerate(df[0]):
        row = proteindata[x]
        data = [row.accession,df[1].iloc[i],row.organism,row.taxonomy,df[2].iloc[i],row.GOTerms,row.sequence]
        csvwriter.writerow(data)
