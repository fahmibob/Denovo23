import requests
import xml.etree.ElementTree as ET

# File containing the list of species, one species per line
species_file = "species_list.txt"
gene = "rbcL"

# Your API key
api_key = "3427cb5832a349c2516ada0aafcccb3b4009"

# Read the list of species from the file
with open(species_file, "r") as f:
    species_list = f.readlines()

# Open a file to save the sequences
with open("sequences.txt", "w") as f:
    for species in species_list:
        species = species.strip() # remove leading and trailing whitespaces

        # Search for gene in species and retrieve specific gene record
        search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term={species}[Organism] AND {gene}[All Fields]&retmax=1&api_key={api_key}"
        search_result = requests.get(search_url).text

        # Extract the ID of the first search result
        root = ET.fromstring(search_result)
        id_list = root.findall("IdList/Id")
        if len(id_list) == 0:
            print(f"No sequence found for {species} and {gene}")
            continue
        id = id_list[0].text

        # Retrieve the sequence of the target gene
        fetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id={id}&rettype=gene_fasta&api_key={api_key}"
        sequence = requests.get(fetch_url).text
        f.write(sequence)
