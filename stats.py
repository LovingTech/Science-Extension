# This is an example of a row is the data.csv file produced the data_extraction.py script.
# The column are as follows: (verbatim)
# - name
# - acession
# - organism
# - taxonomy 
# - sequence
# - db-reference
# -------------------------------------------------------------
# 001R_FRG3G,
# ['Q6GZX4'],
# Frog virus 3 (isolate Goorha) (FV-3).,
# "['Viruses', 'Varidnaviria', 'Bamfordvirae', 'Nucleocytoviricota', 'Megaviricetes', 'Pimascovirales', 'Iridoviridae', 'Alphairidovirinae', 'Ranavirus', 'Frog virus 3']",
# MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVECPKAPVEWNNPPSEKGLIVGHFSGIKYKGEKAQASEVDVNKMCCWVSKFKDAMRRYQGIQTCKIPGKVLSDLDAKIKAYNLTVEGVEGFVRYSRVTKQHVAAFLKELRHSKQYENVNLIHYILTDKRVDIQHLEKDLVKDFKALVESAHRMRQGHMINVKYILYQLLKKHGHGPDGPDILTVKTGSKGVLYDDSFRKIYTDLGWKFTPL,
# ['GO:0046782']
# -------------------------------------------------------------

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import proteintools
import multiprocessing as mp
import seaborn as sns
import scipy.stats as ss
import os
import pickle
import csv


# directory for the graphs

#Initialise the dataframe
df = pd.read_csv('data.csv')

# Results output from PredtoResult.py
results = pd.read_csv('./results/deepgoplus/results_deepgoplus.csv', index_col=0,header=None, names=["Accession", "Similarity","Predicted GO"])
#results1 = pd.read_csv('./results/SPROF/results_sprofFinal.csv', index_col=0,header=None, names=["Accession", "Similarity","Predicted GO"])
#results = pd.read_csv('./results/netgo3/results_netgoFinal.csv', index_col=0, header=None, names=["Accession","Similarity","Predicted GO"])
prefix = "deepgo"

# Organism Frequency Histogram
def organism_barh(num: int = 20):
    plt.barh(df["organism"].value_counts()[:num].index, df["organism"].value_counts()[:num].values)
    plt.show()

def aminoacid_sequencehist():
    data = proteintools.ProteinData("data.csv")
    x = [x.sequence_len for x in data]
    plt.hist(x,bins=5000,density=True,color="black",range=(0,5000))
    plt.xlabel("No of Amino Acids in Sequence")
    plt.ylabel("Frequency")
    plt.savefig(prefix+"_aminoacid_hist.svg")

# Histogram of the amino acid sequence length
def aminoacid_sequence():
    data = proteintools.ProteinData("output.csv")
    x = [data[i].sequence_len for i in results.index]
    y = results[1]
    plt.hist2d(x,y,bins=50,cmap='plasma')
    plt.xlabel("No of Amino Acids in Sequence")
    plt.ylabel("Similarity between Model and SwissProt Dataset")
    plt.savefig(prefix+"_aminoacid_hist2d.svg")

# Produces a violin plot of all similarity across different taxonomies
def taxonomy_violin():
    data = proteintools.ProteinData("output.csv")
    x = [data[i].taxonomy[0] for i in results.index]
    y = results[1]
    taxdf = pd.DataFrame.from_dict({"taxonmy":x,"sim":y})
    violindata = [d["sim"].to_list() for _, d in taxdf.groupby('taxonmy')]
    print([d for _, d in taxdf.groupby('taxonmy')])
    plt.violinplot(violindata, showmeans=True)
    plt.xticks([1,2,3,4],labels=["Archaea","Bacteria","Eukaryota","Viruses"])
    plt.savefig(prefix+"_taxonmy_violins.svg")

# Creates multiple histograms for the different taxonomies and serves as a replacement for `taxonomy_violin`
def taxonomy_hists():
    data = proteintools.ProteinData("output.csv")
    x = [data[i].taxonomy[0] for i in results.index]
    y = results["Similarity"]
    taxdf = pd.DataFrame.from_dict({"taxonomy":x,"sim":y})
    data = [d["sim"].to_list() for _, d in taxdf.groupby('taxonomy')]
    archea = data[0]
    bacteria = data[1]
    eukaryotes = data[2]
    viruses = data[3]
    print(len(archea),len(bacteria),len(eukaryotes),len(viruses))
    fig, axs = plt.subplots(2, 2,sharex=True,sharey=True)
    axs[0, 0].hist(archea,bins=125,color="red", range=(0,1))
    axs[0, 0].set_title("Archaea")
    axs[1, 0].hist(bacteria,bins=125, color="black", range=(0,1))
    axs[1, 0].set_title("Bacteria")
    axs[0, 1].hist(eukaryotes,bins=125, color="green", range=(0,1))
    axs[0, 1].set_title("Eukaryotes")
    axs[1, 1].hist(viruses,bins=125, color="blue", range=(0,1))
    axs[1, 1].set_title("Viruses")
    fig.tight_layout()
    fig.savefig(prefix+"_taxonomy_histograms.svg")

def taxonomy_hist_seaborn():
    data = proteintools.ProteinData("output.csv")
    x = [data[i].taxonomy[0] for i in results.index]
    y = results["Similarity"]
    z = [data[i].sequence_len for i in results.index]
    taxdf = pd.DataFrame.from_dict({"taxonomy":x,"Similarity":y,"Sequence_Length":z})
    data = [d for _, d in taxdf.groupby('taxonomy')]
    archaea = data[0]
    bacteria = data[1]
    eukaryotes = data[2]
    viruses = data[3]
    a = sns.JointGrid(data=archaea, y="Similarity", x="Sequence_Length",ylim=(0,1), space=0)
    a.plot_joint(sns.histplot,color="red")
    a.ax_joint.set(xscale="log")#,log_scale=(True,False))
    a.plot_marginals(sns.histplot,bins=125, color="red")
    a.savefig(prefix+"archaea.svg")

    b = sns.JointGrid(data=bacteria, y="Similarity", x="Sequence_Length",ylim=(0,1), space=0)
    b.plot_joint(sns.histplot, color="black")
    b.ax_joint.set(xscale="log")#,log_scale=(True,False))
    b.plot_marginals(sns.histplot,bins=125, color="black")
    b.savefig(prefix+"bacteria.svg")
    
    e = sns.JointGrid(data=eukaryotes, y="Similarity", x="Sequence_Length",ylim=(0,1), space=0)
    e.plot_joint(sns.histplot, color="limegreen")
    e.ax_joint.set(xscale="log")#,log_scale=(True,False))
    e.plot_marginals(sns.histplot,bins=125,color="limegreen")
    e.savefig(prefix+"eukaryotes.svg")

    v = sns.JointGrid(data=viruses, y="Similarity", x="Sequence_Length",ylim=(0,1), space=0)
    v.plot_joint(sns.histplot,color="royalblue")
    v.ax_joint.set(xscale="log")#,log_scale=(True,False))
    v.plot_marginals(sns.histplot,bins=125,color="royalblue")
    v.savefig(prefix+"viruses.svg")



def add_additional_data(df):
    if os.path.isfile("cache.pickle"):
        with open("cache.pickle", "rb") as file:
            df = pickle.load(file)
    else:
        data = proteintools.ProteinData("output.csv")
        df["taxonomy"] = [data[i].taxonomy[0] for i in df.index]
        df["sequencelen"] = [data[i].sequence_len for i in df.index]
        df["organsim"] = [data[i].organism for  i in df.index]
        df["Labels"] = [data[i].GOTerms for  i in df.index]
        with open("cache.pickle", "wb") as file:
            pickle.dump(df,file)
        df.to_excel(prefix+".xlsx")
    return df

# GO Anotation Frequency Histogram
def no_of_GOA(no_bins: int = 100):
    plt.hist(df["db-references"].apply(lambda x: len(x)),bins=no_bins)
    plt.title("GO Annotation Frequency")
    plt.show()

# Basic Analysis for series of numerical data
def AnalysisForSeries(data):
    print(data.describe())
    print("skew:",data.skew())
    print("kurtosis:", data.kurtosis())
    plt.hist(data,bins=100)
    plt.savefig(prefix+"_serieshist.png")

# number of go anotations distribution. 
# Work in progress ...
def GOA_density_dist():
    df["db-references"].apply(lambda x: len(x))

def samples(df,column_a, value_a, column_b, value_b):
    a_index = df.index[df[column_a] == value_a].tolist()
    b_index = df.index[df[column_b] == value_b].tolist()
    return a_index,b_index

def ttest_viruses_and_eukaryota(results):
    print("---------------") 
    print("viruses and eukaryotes")
    results = add_additional_data(results) 
    x_index, y_index = samples(df=results, column_a="taxonomy", value_a="Viruses", column_b="taxonomy", value_b="Eukaryota")
    x = [results.loc[i]["Similarity"] for i in x_index]
    y = [results.loc[i]["Similarity"] for i in y_index]
    print("Variance x:",np.var(x))
    print("Variance y:",np.var(y))
    #with open("ttesting.csv", "w") as csvfile:
    #    rows = zip(x,y)
    #    writer = csv.writer(csvfile)
    #    for row in rows:
    #        writer.writerow(row)
    print(ss.ttest_ind(x,y,equal_var=True))

def ttest_viruses_and_bacteria(results):
    print("---------------") 
    print("viruses and bacteria")
    results = add_additional_data(results) 
    x_index, y_index = samples(df=results, column_a="taxonomy", value_a="Viruses", column_b="taxonomy", value_b="Bacteria")
    x = [results.loc[i]["Similarity"] for i in x_index]
    y = [results.loc[i]["Similarity"] for i in y_index]
    print("Variance x:",np.var(x))
    print("Variance y:",np.var(y))
    print(ss.ttest_ind(x,y,equal_var=True))

def ttest_viruses_and_archaea(results):
    print("---------------") 
    print("viruses and archaea")
    results = add_additional_data(results) 
    x_index, y_index = samples(df=results, column_a="taxonomy", value_a="Viruses", column_b="taxonomy", value_b="Archaea")
    x = [results.loc[i]["Similarity"] for i in x_index]
    y = [results.loc[i]["Similarity"] for i in y_index]
    print("Variance x:",np.var(x))
    print("Variance y:",np.var(y))
    print(ss.ttest_ind(x,y,equal_var=True))

def ttest_bacteria_and_eukaryota(results):
    print("---------------") 
    print("bacteria and eukaryota")
    results = add_additional_data(results) 
    x_index, y_index = samples(df=results, column_a="taxonomy", value_a="Eukaryota", column_b="taxonomy", value_b="Bacteria")
    x = [results.loc[i]["Similarity"] for i in x_index]
    y = [results.loc[i]["Similarity"] for i in y_index]
    print("Variance x:",np.var(x))
    print("Variance y:",np.var(y))
    print(ss.ttest_ind(x,y,equal_var=True))

def ttest_bacteria_and_archaea(results):
    print("---------------") 
    print("bacteria and archaea")
    results = add_additional_data(results) 
    x_index, y_index = samples(df=results, column_a="taxonomy", value_a="Archaea", column_b="taxonomy", value_b="Bacteria")
    x = [results.loc[i]["Similarity"] for i in x_index]
    y = [results.loc[i]["Similarity"] for i in y_index]
    print("Variance x:",np.var(x))
    print("Variance y:",np.var(y))
    print(ss.ttest_ind(x,y,equal_var=True))

def ttest_eukaryota_and_archaea(results):
    print("---------------") 
    print("eukaryota and archaea")
    results = add_additional_data(results) 
    x_index, y_index = samples(df=results, column_a="taxonomy", value_a="Archaea", column_b="taxonomy", value_b="Eukaryota")
    x = [results.loc[i]["Similarity"] for i in x_index]
    y = [results.loc[i]["Similarity"] for i in y_index]
    print("Variance x:",np.var(x))
    print("Variance y:",np.var(y))
    print(ss.ttest_ind(x,y,equal_var=True))


def precision_recall():
    proteindata = proteintools.ProteinData("output.csv")
    accession_pred = []
    for x in range(len(results)):
        data = results.iloc[x]
        accession_pred.append((data.name,data["Predicted GO"],proteindata[data.name].GOTerms))
    


def modeprediction():
    proteindata = proteintools.ProteinData("output.csv")
    accession_pred = []
    for x in range(len(results)):
        data = results.iloc[x]
        accession_pred.append((data.name,data["Predicted GO"],proteindata[data.name].GOTerms))
    
    true_prediction_count = {}

    for entry in accession_pred:
        _,prediction,groundtruth=entry
        print(groundtruth)
        intersection = list(set(prediction) & set(groundtruth))
        for x in intersection:
            true_prediction_count[x] += 1
    
    return true_prediction_count


# Creating Graphs
def main():
    print("main")
    #ttest_viruses_and_eukaryota(results)
    #ttest_viruses_and_bacteria(results)
    #ttest_viruses_and_archaea(results)
    
    #ttest_bacteria_and_eukaryota(results)
    #ttest_bacteria_and_archaea(results)

    #ttest_eukaryota_and_archaea(results)
    
    #precision_recall()
    #print(modeprediction())

    #print(ss.ttest_ind(results["Similarity"],results1["Similarity"],equal_var=False))
    #print(ss.mannwhitneyu(results["Similarity"],results1["Similarity"]))
    #print(ss.wilcoxon(results["Similarity"][:9000],results1["Similarity"][:9000]))
    #print(ss.bws_test(results["Similarity"],results1["Similarity"]))
    #print(ss.epps_singleton_2samp(results["Similarity"],results1["Similarity"])) # Comparison of underlying distribution
    #print(ss.mood(results["Similarity"],results1["Similarity"])) # Comparison of scale parameters
    #print(ss.kruskal(results["Similarity"],results1["Similarity"])) # Comparison of medians

    
    #taxonomy_hist_seaborn()
    #aminoacid_sequencehist()
    taxonomy_hists()

if __name__ == "__main__":
    main()
