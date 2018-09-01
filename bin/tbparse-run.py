#! /usr/bin/env python

import tbparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--mutation",required=True,help="the name of the database (default='CRyPTIC2')")
    options = parser.parse_args()

    
