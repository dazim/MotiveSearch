#!/usr/bin/env python

################################################################################
##### All imports ##############################################################
################################################################################

import argparse                 # for parsing filename over cli
import itertools                # for generation of n-size k-mers
import matplotlib.pyplot as plt # for visualisation
import pandas as pd             # for kmer count table
import re                       # for regex
from datetime import datetime   # for benchmarking

################################################################################

################################################################################
##### Initial parameters #######################################################
################################################################################

input_type = "default"

################################################################################

################################################################################
##### Filename can be passed to script over the CLI ############################
################################################################################


parser = argparse.ArgumentParser()

parser.add_argument("--mode", "-m", type=str, required=True)
parser.add_argument("--input_genes", "-i", type=str, required=True)
parser.add_argument("--input_type", "-t", type=str, required=False)
parser.add_argument("--genome", "-g", type=str, required=False)

args = parser.parse_args()


################################################################################

################################################################################
##### Extracts all C.gl. identifiers from any input ############################
################################################################################

def extract_counts(inputfile, inputtype):
    file = open(inputfile, "r")

    gene_identifiers = []

    if inputtype == "counttable":

        for line in file:

            try:

                row = re.findall("cg\d{4}.+", line)[0]
                parts = str(row).split("\t")
                gene_identifiers.append([parts[0], parts[2], parts[3], parts[4],
                                         parts[6]])

            except IndexError:

                pass  # print("Error in extract_genes")

    else:

        for line in file:
            gene_identifiers.append(re.findall("cg\d{4}", str(line)))

    return gene_identifiers


################################################################################

################################################################################
##### Load promoter sequences ##################################################
################################################################################

def load_promoters(inputfile):
    file = open(inputfile, "r")

    promoters = []

    for line in file:

        try:

            line = str(line).strip().split("\t")
            # parts = str(row).split("\t")
            if line[4] == "fwd":
                strand = "+"
            else:
                strand = "-"

            promoters.append([line[0], line[2], line[1], strand, line[5]])

        except IndexError:

            print ("Error in load_promoters")

    return promoters


################################################################################


################################################################################
##### Extract n most frequent kmers from used promoters ########################
################################################################################

def top_kmers(promoters, kmers, n):
    seq = ""

    df = pd.DataFrame(index=kmers, columns=["Total"])

    for promoter in promoters:
        seq += promoter[4]

    for i in range(0, len(kmers)):
        df.at[kmers[i], "Total"] = seq.count(kmers[i])

    sorted_df = best_kmers(df)

    return sorted_df.index.tolist()[0:n]


################################################################################

################################################################################
##### Extract promoters of DEGs ################################################
################################################################################

def count_kmers(promoters, kmers):
    gen_identifier = []

    for promoter in promoters:
        gen_identifier.append(promoter[1])

    gen_identifier.append("Total")

    df = pd.DataFrame(index=kmers, columns=[gen_identifier])

    for i in range(0, len(promoters)):

        for j in range(0, len(kmers)):
            df.at[kmers[j], promoters[i][1]] = promoters[i][4].count(kmers[j])

    for m in range(0, len(kmers)):

        for n in range(0, len(promoters)):
            df.at[kmers[m], "Total"] = 0

    for w in range(0, len(kmers)):

        for v in range(0, len(promoters)):
            df.at[kmers[w], "Total"] += df.at[kmers[w], promoters[v][1]]

    return df


################################################################################

################################################################################
##### Count total appearences of kmers #########################################
################################################################################

def summarize_kmers(df):
    summary = []

    return summary


################################################################################

################################################################################
##### Extract promoters of DEGs ################################################
################################################################################

def filter_promoters(promoters, genes):
    filtered_promoters = []

    for gene in genes:

        for i in range(0, len(promoters)):

            if gene[0] == promoters[i][1]:
                filtered_promoters.append(promoters[i])

    return filtered_promoters


################################################################################

################################################################################
##### Generate kmers of specified length #######################################
################################################################################

def generate_kmers(size):
    bases = ["a", "g", "t", "c"]

    kmers = [''.join(i) for i in itertools.product(bases, repeat=size)]

    return kmers


################################################################################

################################################################################
##### Cout frequency of kmers over different promoters #########################
################################################################################

def best_kmers(df):
    return df.sort_values("Total", ascending=0)


################################################################################

################################################################################
##### Cout frequency of kmers over different promoters #########################
################################################################################

def visualize_kmers(df):
    plt.bar(range(len(df["Total"].values)), df["Total"].values)

    plt.show()

    return


################################################################################

################################################################################
##### Cout frequency of kmers over different promoters #########################
################################################################################

def prepare_genome(file, kmers):
    file = open(file, "r")

    gene_identifiers = []

    seq = ""

    df = pd.DataFrame(index=kmers, columns=["Count"])

    for line in file:
        seq += line.strip().lower()

    for i in range(0, len(kmers)):
        df.at[kmers[i], "Count"] = seq.count(kmers[i])

    return df


################################################################################



################################################################################
##### Sandbox ##################################################################
################################################################################


if args.mode == "AnalysePromoters":

    gene_counts = extract_counts(args.input_genes, args.input_type)
    print(len(gene_counts))

    dge_genes = [["cg0002"], ["cg0012"], ["cg0003"], ["cg0004"], ["cg0005"]]

    promoters = load_promoters("promoters_cg_short.txt")
    print(len(promoters))

    c = filter_promoters(promoters, dge_genes)

    d = generate_kmers(3)

    e = count_kmers(c, d)

    # print(e)

    # f = visualize_kmers(e)

    f = best_kmers(e)

    visualize_kmers(f)

elif args.mode == "PrepareGenome":

    t1 = datetime.now()

    gene_counts = extract_counts(args.input_genes, args.input_type)
    print(len(gene_counts))

    dge_genes = [["cg0002"], ["cg0012"], ["cg0003"], ["cg0004"], ["cg0005"]]

    promoters = load_promoters("promoters_cg_short.txt")
    print(len(promoters))

    c = filter_promoters(promoters, dge_genes)

    kmers_best = top_kmers(c, generate_kmers(10), 20)

    for kmer in kmers_best:
        print kmer

#    for i in range(0, 21):
#        t1 = datetime.now()

        # print(prepare_genome(args.input_genes, kmers_best))


#        print(str(i), " ", str(datetime.now() - t1))
# print(summarize_kmers(e))



################################################################################
