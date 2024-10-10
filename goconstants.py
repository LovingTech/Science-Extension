import neo4j

def getconstants(URL, AUTH):
    driver = neo4j.GraphDatabase.driver(uri=URL, auth=AUTH)
    # All molecular function desecendents calculated
    global MF_Descentents 
    MF_Descentents= driver.execute_query(
            """MATCH p=(:Class{name: "GO_0003674"})<-[:is_a*]-(:Class)
            UNWIND nodes(p) AS nodes
            RETURN count(DISTINCT nodes) AS result""",
            result_transformer_=neo4j.Result.to_df
            )["result"][0]
    # All biological process desecendents calculated
    global BP_Descentents
    BP_Descentents = driver.execute_query(
            """MATCH p=(:Class{name: "GO_0008150"})<-[:is_a*]-(:Class) 
            UNWIND nodes(p) AS nodes
            RETURN count(DISTINCT nodes) AS result""",
            result_transformer_=neo4j.Result.to_df
            )["result"][0]
    # All cellular component desecendents calculated
    global CC_Descentents 
    CC_Descentents = driver.execute_query(
            """MATCH p=(:Class{name: "GO_0005575"})<-[:is_a*]-(:Class) 
            UNWIND nodes(p) AS nodes 
            RETURN count(DISTINCT nodes) AS result""",
            result_transformer_=neo4j.Result.to_df
            )["result"][0]
    
    global Total_Terms
    # Commented out as I considered it inconsistent
    # 
    #Total_Terms = driver.execute_query(
    #        """MATCH (n:Class) RETURN count(*) AS result""", 
    #        result_transformer_=neo4j.Result.to_df
    #        )["result"][0]
    
    Total_Terms = CC_Descentents + BP_Descentents + MF_Descentents
    driver.close()

# TESTING
if __name__ == "__main__":
    import personal
    getconstants(personal.URI, personal.AUTH)
    print("MF_Descentents:", MF_Descentents)
    print("BP_Descentents:", BP_Descentents)
    print("CC_Descentents:", CC_Descentents)
    print("Total_Terms:", Total_Terms)
