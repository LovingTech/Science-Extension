# Similarity measures derived from
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-S5-S4

import GOToolsNeo4j

class TermSemanticSimilarity():
    def __init__(self, TermA: GOToolsNeo4j.GOTerm, TermB: GOToolsNeo4j.GOTerm):
        # Create a method to print the Terms
        self.TermA = TermA
        self.TermB = TermB
    
    ## Resnik Similarity
    def Resnik(self):
        # Get the IC of the most informative common ancestor (TermB)
        Ca , _, _ = GOToolsNeo4j.commenancestor(self.TermA, self.TermB)
        # Get the most informative common ancestor
        # Get the IC of the most informative common ancestor
        # Calculate the Resnik Similarity
        if Ca is None:
            return None
        return Ca.InformationContent()
    
    ## Lin's Term semantic similarity
    def Lin(self):
        # Get the common ancestor
        Ca , _, _ = GOToolsNeo4j.commenancestor(self.TermA, self.TermB)
        if Ca is None:
            return None
        IC_Ca = Ca.InformationContent()
        # Get the IC for both terms
        IC_Ta = self.TermA.InformationContent()
        IC_Tb = self.TermB.InformationContent()
        # The formula for Lin's measure
        similarity = (2 * IC_Ca)/(IC_Ta + IC_Tb)
        return similarity
    
    def JiangAndConrath(self):
        # Get the common ancestor
        Ca , _, _ = GOToolsNeo4j.commenancestor(self.TermA, self.TermB)
        if Ca is None:
            return None
        IC_Ca = Ca.InformationContent()
        # Get the IC for both terms
        IC_Ta = self.TermA.InformationContent()
        IC_Tb = self.TermB.InformationContent()
        # The formula for Lin's measure
        similarity = 1 + IC_Ca - ((IC_Ta + IC_Tb)/2)
        return similarity

class ProteinFunctionSimilarity:
    def __init__(self,ProteinPred:list[GOToolsNeo4j.GOTerm],ProteinLabel:list[GOToolsNeo4j.GOTerm]) -> None:
        self.ProteinPred = ProteinPred
        self.ProteinLabel = ProteinLabel
    def Average(self):
        # Get the average semantic similarity between two proteins
        similarity = 0
        terms_compared = len(self.ProteinPred)*len(self.ProteinLabel)
        for TermA in self.ProteinPred:
            for TermB in self.ProteinLabel:
                 simterm = TermSemanticSimilarity(TermA,TermB).Lin()
                 print(TermA,TermB,simterm)
                 if simterm is None:
                     terms_compared -= 1
                 else:
                     similarity += simterm
        if terms_compared == 0:
            return 0        
        else:
            return similarity/terms_compared

    def Max(self):
        # Get the max semantic similarity between two proteins
        score = []
        for TermA in self.ProteinPred:
            for TermB in self.ProteinLabel:
                simterm = TermSemanticSimilarity(TermA,TermB).Lin()
                score.append(simterm)
        return max(score)

        

# TESTING
if __name__ == "__main__":
    # Create a GOTerm object
    known = list(map(GOToolsNeo4j.GOTerm,['GO:0016020', 'GO:0110165', 'GO:0005575']))
    prediction = list(map(GOToolsNeo4j.GOTerm,['GO:0110165', 'GO:0005575', 'GO:0008150', 'GO:0003674', 'GO:0005622', 'GO:0016020', 'GO:0005737', 'GO:0009987', 'GO:0008152', 'GO:0071944', 'GO:0044237', 'GO:0071704', 'GO:0003824']))
    print(ProteinFunctionSimilarity(prediction,known).Average())
    #TermA = GOToolsNeo4j.GOTerm("GO:0005542")
    #TermB = GOToolsNeo4j.GOTerm("GO:0005542")
    #TermB = GOToolsNeo4j.GOTerm("GO:0030170")
    # Create a SemanticSimilarity object
    #Sim = TermSemanticSimilarity(TermA, TermB)
    #conancestor, _, _ = GOToolsNeo4j.commenancestor(TermA, TermB,visual=True)
    # Calculate the Resnik Similarity
    #print("Resnik: " + str(Sim.Resnik()))
    #print("Lin: " + str(Sim.Lin()))
    #print("Jiang and Conrath: " + str(Sim.JiangAndConrath()))
    
