
# Running the Code
`requirements.txt` will have the required libraries and versions. Which can be installed using `pip install -r requirements.txt`

## data_extraction
If this file is run as main then it will create a file with the name `testing.csv` using the `get_data()` function.
### `get_data() -> None`
**Notes**
This function will first search for the file given in the `filepath` parameter and if it cannot find the file or parse it correctly the function will automatically attempt to download the `uniprot_sprot.data` file, unpack it and parse that. Then it will write data into a csv file with path in `outfile`. With the rows `name,accession,organism,taxonomy,sequence,db-references`. The last of which will only contain Gene Ontology references.
**Parameters**
- `filepath: str`
	- This is used to show where the `uniprot_sprot.dat` file is. Using either absolute or relative path.
	- Default: `uniprot_sprot.data`
- `outfile: str`
	- Used when writing the file that contains all the relevant data.
## protein-selection
It assumes that the protein data is stored in a file called `data.csv` and will other wise call functions from `data_extraction` to get that data, going as far as to re-download the entire `uniprot_sprot.dat` file if it cannot find it. 

### Running as a script
When being run as a script it has parameters that can be seen by using `python protein-selection.py --help`. This should provide some idea of how to use it.

### `parse_string_to_list() -> list[str]`
**Notes**
Basically it takes a string of a list and turns it into list, this function isn’t perfect because if one of strings in the list contain a comma it will not parse it correctly, but considering what is being parse this is not a problems and _should not be used for anything_ but inside of already existing code.
### Class `Row()`
**Notes**
This class was created to make the code significantly more readable and perform functions in a more concise manner. 
**Parameters**
* `no: int`
	* Integer of row index in `proteindata` dataframe

### `get_proteins() -> list`
**Notes**
Collects a specific subset of proteins that fits a specific selection criteria. The logic behind this function isn’t optimal as it will continue to select at random ones that have already be selected and/or those that aren’t eligible before check that they aren’t. 
**Parameters**
* `proteindata` 
	* is `pandas.dataframe` object that contains all of the required data
* `number: int`
	* Used to tell how many proteins should be collected 
* `GOterm_Label_No: int`
	* The minimum number of Gene ontology terms a protein can have to be selected for comparison
* `csv_out: str`
	* The name of the file that the selected proteins will be contained in. 
### `validation()`
**Notes**
This was used in testing to confirm that the data from the file is consistent with what was meant to be selected. Should not be required for normal use.
**Parameters**
* `path: str`
	* Used for the path to the file to check the values of.

## proteintools
**Notes**
Used in the comparison to quickly get relevant data about each protein
