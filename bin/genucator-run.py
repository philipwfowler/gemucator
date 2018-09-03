#! /usr/bin/env python

import argparse
from genucator import genucator

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--mutation",required=True,help="the name of the database (default='CRyPTIC2')")
    options = parser.parse_args()

    tb_reference_genome=genucator()

    (locations,bases)=tb_reference_genome.locate_mutation(options.mutation)

    print(options.mutation+":")
    for (p,b) in zip(locations,bases):
        print(p,b)
