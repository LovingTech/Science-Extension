import requests

# Grabs graphs of the ontology from QuickGO
def GetGraph(ids:list|str, ShowChildren:bool=False, ShowKey:bool=True):
    # Used for the api request for the chart
    api_base = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{ids}/chart?ids={ids}&showKey={showkey}&showChildren={showchildren}"
    # Creates a title for the image
    title = f"./graphs/Graph of {" and ".join(ids).replace("/","_")}.png"
    # Allows the function to deal with multiple IDs in one graph
    if isinstance(ids,list):
        ids = ",".join(ids)
    url = api_base.format(ids=ids, showkey=ShowKey, showchildren=ShowChildren)
    response = requests.get(url,headers={ "Accept" : "image/png"})
    if not response.ok:
        response.raise_for_status()
    with open(title, "wb") as f:
        # Write the data for the image
        f.write(response.content)

# TESTING
if __name__ == "__main__":
    #GetGraph(["GO:0004386", "GO:0140490", "GO:0033856"], ShowChildren=True, ShowKey=True)
    #GetGraph(["GO:0043231", "GO:0005488"], ShowChildren=True, ShowKey=True)
    GetGraph(["GO:0070288", "GO:0140691","GO:0015979"], ShowChildren=True, ShowKey=True)
