from Bio import Entrez, Medline


def search_terms(term):
    """
    Searching the PubMed with the query.
    :param term: The given query.
    :return: Count: for the amount of queries that has been found.
    Id_list for the PudMed ID's.
    """
    Entrez.email = "cca.vanharen@student.han.nl"
    handle = Entrez.esearch(db="pubmed", term=term, retmax=100)
    record = Entrez.read(handle)
    handle.close()
    count = (record["Count"])
    id_list = (record["IdList"])

    return count, id_list


def making_query(compounds, genes, molecular):
    """
    Combining the 3 files into queries.
    :param compounds: The compounds file
    :param genes: The genes file
    :param molecular: The molecular_effects file
    :return: id_list_most could be used in the next function
    """
    for compound in compounds.split("\n"):
        # print(compound)
        for gen in genes.split("\n"):
            # print(gen)
            for molecule in molecular.split("\n"):
                # print(molecule)

                query = "Query: (" + compound + " [tiab] AND " + gen + \
                        " [tiab] AND " + molecule + " [tiab])"
                print(query)
                count, id_list = search_terms(query)
                most = 0
                id_list_most = []
                if int(count) > most:
                    id_list_most = id_list
                get_article(id_list_most)


def get_article(id_list):
    """
    With the known PudMed ID's searching for the title, abstract and
    the author of the article.
    :param id_list: The list with the PubMed ID's
    :return: The info list with the info of the article if there is a match.
    """
    for pubmed_id in id_list:
        search_handle = Entrez.efetch(db="pubmed",
                                      id=pubmed_id,
                                      rettype="medline",
                                      retmode="text",
                                      sort="first+author")
        records = list(Medline.parse(search_handle))
        title = records[0].get("TI", "?")
        abstract = records[0].get("AB", "?")
        author = records[0].get("FAU", "?")
        if title != "?" or abstract != "?" or author != "?":
            info = [pubmed_id, title, abstract, author]
            print(info)
        else:
            print("There is an item missing.")


def main():
    # "Cypripedioideae [ab] AND matK [ab]"
    compounds = open("small_compounds.txt").read()
    genes = open("small_genes.txt").read()
    molecular = open("small_molecular_effects.txt").read()
    making_query(compounds, genes, molecular)
    count, id_list = search_terms("Cypripedioideae [tiab] AND matK ["
                                  "tiab]")
    get_article(id_list)


main()
