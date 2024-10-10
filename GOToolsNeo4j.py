import json
import neo4j
import numpy as np
import requests
import personal
import goconstants
import GOOntologyTools

URI = personal.URI
AUTH = personal.AUTH

driver = neo4j.GraphDatabase.driver(uri=URI, auth=AUTH)
#Get constants
goconstants.getconstants(URI, AUTH)

# Create GOTerm Object
class GOTerm:
    def __init__(self, ID):
        self.ID = ID
        # In the database they replace the colon with a underscore
        self.dbID = ID.replace(":", "_")
        # Creates link directly to the Quickgo website to verify information and get more.
        self.QuickGOUrl = "https://www.ebi.ac.uk/QuickGO/term/" + self.ID

    # Create a method to print the Term:s
    def __str__(self) -> str:
        return self.ID
    # Checks whether a term is a child of self
    def __contains__(self,term):
        if hasattr(term, 'dbID') is True:
            containquery = neo4j.Query("""
                                        MATCH (n:Class{name: $term })
                                        MATCH (m:Class{name: $self })
                                        RETURN exists((n)-[:is_a*]->(m)) AS result
                                       """)
            result = dbquery(containquery,{"term":term.dbID, "self":self.dbID})["result"][0] 
        else:
            result = False
        return result

    def Descendents(self):
        # Query for descentents of a term
        descquery = neo4j.Query("""
         MATCH p=(:Class{name: $term })<-[*]-(:Class)
         UNWIND nodes(p) AS nodes 
         RETURN count(DISTINCT nodes) AS result
         """)
        # Runs the query
        result = dbquery(descquery, {"term":self.dbID})["result"][0] #but it shouldn't return mutliple
        return result

    def Ancestors(self):
        ancestquery = neo4j.Query("""
        MATCH (term:Class{name: $term})
        MATCH (namespace:Class)
        MATCH p=(term)-[*]->(namespace)
        WHERE namespace.name IN ["GO_0003674","GO_0008150","GO_0005575"]
        with collect(nodes(p)) as listOflistOfnodes
        unwind listOflistOfnodes as list
        unwind list as element
        with collect(distinct element.name) as ancestors
        RETURN ancestors
        """)
        # Runs the query
        result = dbquery(ancestquery, {"term":self.dbID})["ancestors"][0]
        resultf = [sub.replace('_', ':') for sub in result]
        return resultf

    # Returns the information content of the term which is required for similarity calculations. 
    # Slower and not optimised
    def InformationContent_WithRespectToDepth(self) -> float:
        icquery = neo4j.Query("""
                MATCH (term:Class{name: $term })
                MATCH (namespace:Class)
                WHERE EXISTS {
                    MATCH (term)-[*]->(namespace)
                    WHERE namespace.name IN ["GO_0003674","GO_0008150","GO_0005575"]
                }
                MATCH p=(term)-[*]->(namespace)
                RETURN namespace.name AS namespace, length(p) as depth,COUNT {
                    MATCH d=(:Class)-[:is_a*]->(term)
                    UNWIND nodes(d) as nodes
                    RETURN nodes} AS descendants
                ORDER BY length(p)
                LIMIT 1
                """)
        result = dbquery(icquery,{"term":self.dbID}).iloc[0]
        # Calculates the information of the term using a topology-based method
        IC = (1 - (np.log(result["descendants"] + 1)/np.log(goconstants.Total_Terms)))*result["depth"]
        return IC
    
    # Get the information content simply
    def InformationContent(self) -> float:
        desc = self.Descendents()
        IC = (1 - (np.log(desc + 1)/np.log(goconstants.Total_Terms)))
        return IC
    
    # Gets any extra data about the GO terms from QuickGO API
    def ExtraInfo(self, local: bool = False) -> dict:
        datadict = {}
        if local is False:
            print("Using the Online Method. Beware when using this in for loops as it can be extremely slow")
            # API Call to get all remaining information about GOTerm, will be useful for when
            # reviewing cousin terms and their common ancestor
            RequestUrl = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{
                self.ID}/complete"
            r = requests.get(RequestUrl, headers={
                             "Accept": "application/json"})
            # Check that it was successful if not then fail safely
            if not r.ok:
                r.raise_for_status()
            # Parse the information
            responseBody = json.loads(r.text)["results"][0]
            # Store in a Python readable format, dictionaries are ideal for this purpose
            datadict = {
                "name": responseBody["name"],
                "definition": responseBody["definition"]["text"],
                "aspect": responseBody["aspect"],
                "isObsolete": responseBody["isObsolete"],
            }
            # Stores as a local variable for later use is nesscsary, therefore removing any repeats
            # of unecesary work.
            self.ExtraData = datadict
        return datadict


def IC_list(goterms):
    query = neo4j.Query("""
                        WITH $list AS terms
                        MATCH (p)
                        WHERE p.name IN terms
                        MATCH (n:Class)-[:is_a*]->(p)
                        WITH COUNT (DISTINCT n) as descendants,terms
                        WITH COLLECT(1-(log(descendants+1))/(log(43121))) as ICp
                        RETURN ICp
                        """)
    result = dbquery(query,{"list":goterms.replace(":","_")})["ICp"][0]
    return result


# Basic Wrapper for the Neo4j Query to ensure that the database isn't accidently 
# written to and returns the results in a nice way.
def dbquery(query: neo4j.Query, params:dict):
    result = driver.execute_query(
        query, 
        parameters_=params,
        routing_="r", 
        database_="neo4j",
        result_transformer_=neo4j.Result.to_df
    )
    return result

# Takes a list of GO term ids and turn them into a list of GO Term objects
def list_to_terms(iter:list) -> list[GOTerm]:
    return list(map(GOTerm,iter))
    
# Get the closest common ancestor between the two terms
def commenancestor(TermA: GOTerm, TermB: GOTerm,visual:bool = False):
    # Checks if the terms are the same
    if TermA.ID == TermB.ID:
        return TermA, 0, 0
    # Creates the query for the common ancestor of both terms and tells the distance 
    # each term is from the common ancestor
    comancestquery = neo4j.Query("""
    MATCH (c1:Class{name: $terma })
    MATCH (c2:Class{name: $termb })
    MATCH (ancestor:Class)
    WHERE EXISTS{
        MATCH (c1)-[*]->(ancestor)
        MATCH (c2)-[*]->(p2:Class)
        WHERE ancestor.name = p2.name
    }
    MATCH pathA = (c1)-[:is_a*]->(ancestor)
    MATCH pathB = (c2)-[:is_a*]->(ancestor)
    RETURN ancestor, length(pathA), length(pathB)
    ORDER by length(pathA)+length(pathB)
    LIMIT 1
    """)
    # Runs the query
    try:
        ancestordbid, lenA, lenB = dbquery(comancestquery,{"terma": TermA.dbID,"termb": TermB.dbID}).to_records(index=False)[0] # Should only be one result because there is a LIMIT 1 in the cypher query
    except: 
        # deals with terms that don't have a common ancestor
        return None, None, None 
    # Returns the ancestor as a GOTerm object
    ancestor = GOTerm(ancestordbid["name"].replace("_",":"))
    # Gets a graph from the QuickGO api to allow for easy visualisation, Has been adjusted to give decent looking results
    if visual is True:
        GOOntologyTools.GetGraph([ancestor.ID,TermA.ID,TermB.ID])
    return ancestor, lenA, lenB


# TESTING
def main():
    #BP = GOTerm("GO:0008150")
    #MF = GOTerm("GO:0003674")
    TERM = GOTerm("GO:0005759")
    #print(TERM.Ancestors())
    #print(TERM.InformationContent())
    #print(BP.InformationContent())
    #print(MF.InformationContent())
    #print(goconstants.Total_Terms)
    #commonancestor, lenA, lenB = commenancestor(TERM, TERM,visual=False)
    #print(commonancestor)
    #print(commonancestor.InformationContent())
    #print(TERM.InformationContent())
    #print(TERM.Ancestors())


if __name__ == "__main__":
    main()
