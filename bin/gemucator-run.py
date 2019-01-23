#! /usr/bin/env python

import argparse
from gemucator import gemucator

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--mutation",help="the full length mutation to locate on the genome e.g. rpoB_S450L")
    parser.add_argument("--nucleotide",action='store_true',help="overrides the logic and forces the code to treat this mutation as a simple nucleotide mutation. Useful for RNA genes.")
    parser.add_argument("--location",type=int,help="the genome position we want to find the gene for")
    parser.add_argument("--promoter_length",type=int,default=100,help="the number of bases upstream of a gene that we will consider to form the promoter")
    parser.add_argument("--genbank_file",help="the path to the reference GenBank file we want to work with. If not specified, the code will load H37rV.gbk from its config/ folder.")
    options = parser.parse_args()

    if options.genbank_file:
        reference_genome=gemucator(genbank_file=options.genbank_file)
    else:
        reference_genome=gemucator()

    if options.mutation:

        (result,data)=reference_genome.locate_mutation(options.mutation,nucleotide_mutation=options.nucleotide)

        if result is True:
            print(options.mutation+":")
            for (p,b) in zip(data[0],data[1]):
                print(b,p)
        else:
            raise ValueError("Mutation "+options.mutation+" does not validate against the supplied GenBank file because the reference amino acid/base is "+data[0]+" !")

    elif options.location:

        (gene,ref,position) = reference_genome.identify_gene(options.location,promoter_length=options.promoter_length)

        if gene is not None:

            print(gene+"_"+ref+str(position))

        else:

            print("Cannot find a gene or plausible promoter for location "+str(options.location))

    else:
        raise Exception("Must specify one of --mutation or --location!")
