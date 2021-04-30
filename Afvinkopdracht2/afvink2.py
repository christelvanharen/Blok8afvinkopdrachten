import sys
from Bio import Entrez
from Bio import Medline
import matplotlib.pyplot as plt
import numpy as np

def search_terms():
    """
    Vraagt naar hoeveel zoektermen gezocht moeten worden en zoekt
    daarop.
    @return: de ingevoerde zoektermen
    """
    terms = []
    try:
        terms_count = int(input("Amount of searchterms: "))
        for i in range(terms_count):
            term = input("Seacrhterm: ")
            terms.append(term)
    except ValueError:
        print("Not a valid number.")
        sys.exit()
    return terms

def id(term):
    """
    Zoekt in de pubmed naar de opgegeven zoektermen.
    @param term: de meegegeven zoekterm
    @return: lijst met gevonden id's in de pubmed
    """
    Entrez.email = "cca.vanharen@student.han.nl"
    handle = Entrez.esearch(db="pubmed", term=term, retmax="1000")
    record = Entrez.read(handle)
    handle.close()
    idlist = record["IdList"]
    print(record)
    return idlist

def get_years(idlist):
    """
    Zoekt naar de jaartallen bij de gevonden id's
    @param idlist: de gevonden id's uit de pubmed
    @return: de lijst met jaren
    """
    handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline",
                           retmode="text")
    records = Medline.parse(handle)
    records = list(records)
    years = []
    for record in records:
        try:
            years.append(int(record.get("DP").split(" ")[0]))
        except AttributeError:
            pass
        except ValueError:
            years.append(int(record.get("DP").split(" ")[1]))
    return years


def get_start_stop_year(years, pos):
    """
    Haalt de range op van de jaren
    @param years: de gevonden jaartallen bij de id's
    @param pos: de positie van de jaartallen
    @return: de jaartallen
    """
    if years[pos] % 10 > 5:
        year = years[pos] - years[pos] % 10 + 5
    else:
        if pos == -1:
            year = years[pos] - years[pos] % 10 + 5
        elif pos == 0 and years[0] % 10 == 0:
            year = years[pos] - 5
        else:
            year = years[pos] - years[pos] % 10
    # print(year)
    return year


def make_year_count_dict(start_year, stop_year):
    """
    De jaren per 5 jaar bepalen
    @param start_year: eerste jaar in de lijst
    @param stop_year: laatste jaar in de lijst
    @return: dictionary, keys; de jaren per 5, values; 0
    """
    year_count = {}
    for i in range(start_year + 1, stop_year, 5):
        year_count["{}-{}".format(i, i + 4)] = 0
    return year_count


def count_years(years, year_count):
    """
    Telt het aantal jaren dat bij die 5 jaar hoort
    @param years: de gevonden jaartallen bij de id's
    @param year_count: dictionary, keys; de jaren per 5, values; 0
    @return: dictionary, keys; de jaren per 5 jaar, values; het
    aantal jaren
    """
    for year in years:
        for key, value in year_count.items():
            if int(key.split("-")[1]) >= year >= int(key.split("-")[0]):
                year_count[key] += 1
    return year_count


def get_unique_keys(year_counts):
    """
    Haalt de unieke keys op uit de year_counts dictionary.
    @param year_counts: dictionary, keys; de jaren per 5, values; 0
    @return: unieke keys in een set
    """
    year_keys = []
    for year_count in year_counts:
        for key in year_count:
            year_keys.append(key)
    year_key_set = sorted(set(year_keys))
    return year_key_set


def hits_in_list(year_counts, year_key_set):
    """
    Plaatst de gevonden hits in een lijst
    @param year_counts: dictionary, keys; de jaren per 5, values; 0
    @param year_key_set: unieke keys in een set
    @return: de gevonden hits in een lijst
    """
    hits_counts = []
    for year_count in year_counts:
        hits_count = []
        for year_key in year_key_set:
            # Als de key niet in de dictionary zit
            if year_count.get(year_key) is None:
                hits_count.append(0)
            # Als de key wel in de dictionary zit
            else:
                hits_count.append(year_count.get(year_key))
        hits_counts.append(hits_count[:])
        hits_count.clear()
    return hits_counts


def bar_plot(hits_counts, year_key_set, terms_with_results):
    """
    Maakt de grafiek, met het aantal artikelen op de y-as en het
    aantal jaren per 5 op de x-as.
    @param hits_counts: de gevonden hits in een lijst
    @param year_key_set: unieke keys in een set
    @param terms_with_results: termen met resultaten
    @return: staafdiagram
    """
    x_pos = [np.arange(len(hits_counts[0]))]
    for i in range(len(year_key_set)):
        x_pos.append([x + 0.25 for x in x_pos[i]])

    for i in range(len(hits_counts)):
        plt.bar(x_pos[i], hits_counts[i], width=0.25,
                label=terms_with_results[i])

    plt.xlabel("In jaren")
    plt.ylabel("Aantal artikelen")
    plt.title("Aantal artikelen per vijf jaar")
    plt.xticks([r + 0.25 for r in range(len(hits_counts[0]))],
               year_key_set)
    plt.xticks(size=9)
    plt.legend()
    plt.show()


def main():
    terms = search_terms()
    terms_with_results = []
    year_counts = []

    for term in terms:
        idlist = id(term)
        years = get_years(idlist)
        if len(years) == 0:
            print("Geen resultaten gevonden voor: {}".format(term))
        else:
            terms_with_results.append(term)
            years.sort()
            start_year = get_start_stop_year(years, 0)
            stop_year = get_start_stop_year(years, -1)
            year_count = make_year_count_dict(start_year, stop_year)
            year_count = count_years(years, year_count)
            year_counts.append(year_count)

    if len(terms_with_results) == 0:
        sys.exit()
    else:
        year_key_set = get_unique_keys(year_counts)
        hits_counts = hits_in_list(year_counts, year_key_set)
        bar_plot(hits_counts, year_key_set, terms_with_results)


main()
