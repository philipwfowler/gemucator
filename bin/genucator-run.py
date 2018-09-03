#! /usr/bin/env python

import argparse
from genucator import genucator

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--mutation",help="the full length mutation to locate on the genome e.g. rpoB_S450L")
    parser.add_argument("--nucleotide",action='store_true',help="overrides the logic and forces the code to treat this as a simple nucleotide mutation. Useful for RNA genes.")
    parser.add_argument("--location",type=int,help="the genome position we want to find the gene for")
    options = parser.parse_args()

    tb_reference_genome=genucator()

    if options.mutation:

        (locations,bases)=tb_reference_genome.locate_mutation(options.mutation,nucleotide_mutation=options.nucleotide)
            # (locations,bases)=tb_reference_genome.locate_mutation(options.mutation)

        print(options.mutation+":")
        for (p,b) in zip(locations,bases):
            print(p,b)

    elif options.location:

        (gene,ref,position) = tb_reference_genome.identify_gene(options.location)

        print(gene+"_"+ref+str(position))

    else:
        raise Exception("Must specify one of --mutation or --location!")
