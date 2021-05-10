from Bio import Entrez

def search_terms():
    term = ("")
    # "Cypripedioideae [ab] AND matK [ab]",
    Entrez.email = "cca.vanharen@student.han.nl"
    handle = Entrez.esearch(db="pubmed", term=term, retmax=100)

    # Reads the esearch and saves the ID's in id_list
    record = Entrez.read(handle)
    handle.close()
    id_list = (record["IdList"])
    print(id_list)

    return id_list

def terms():


def get_pubtator_link(id_list):
    """This function uses the id's from the id_list to create a pubtator link
    which filters the genes, mutations and diseases out of the title and
    abstract.

    :param id_list: List with all found pubmed id's.
    :return pubtator_link: Pubtator link with the title, abstact, genes,
    diseases and mutations of each article.
    """
    # Standard pubtator link format:
    link = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications" \
           "/export/pubtator?pmids=idvalues&concepts=gene,mutation,disease"

    # Creates a string with all the ID's. Separated by a comma.
    id_string = ""
    for i in range(len(id_list)):
        id_string += id_list[i] + ","
    id_string = id_string[:-1]

    # Adds the ID values to the standard link.
    pubtator_link = link.replace("idvalues", id_string)
    print(pubtator_link)

    return pubtator_link

def main():
    compounds = open("small_compounds.txt")
    genes = open("small_genes.txt")
    molecular = open("small_compounds.txt")
    id_list = search_terms()
    get_pubtator_link(id_list)

main()