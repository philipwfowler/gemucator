[![DOI](https://zenodo.org/badge/147030278.svg)](https://zenodo.org/badge/latestdoi/147030278)

# gemucator

gemucator is short for "Genbank Mutation Locator". It is a simple Python3 class for incorporating into tools that, if you give it a mutation, it will tell its location in the reference genome (and vice versa). The `gemucator` class accepts a path to a genbank file; since I am working with *M. tuberculosis* this is the H37rV genbank file by default, but any genbank file should work.

The package comes with a simple script called `gemucator-run.py` that shows how it works. All these examples are for TB.

First, let's see what happens when we give it an amino acid mutation (which has to be, be definition, in the coding sequence of a gene).

```
> gemucator-run.py --mutation rpoB_S450L
rpoB_S450L:
761153 t
761154 c
761155 g
```
It returns three rows, since there are three bases in the triplet, each with the position in the H37rV reference genome and the reference base. Note that the code is very defensive and checks that `tcg` is a Serine, which happily it is. If we get the mutation wrong, then the code will stop and catch fire.

```
> gemucator-run.py --mutation rpoB_K450L
Traceback (most recent call last):
  File "/Users/fowler/Library/Python/3.5/bin/gemucator-run.py", line 6, in <module>
    exec(compile(open(__file__).read(), __file__, 'exec'))
  File "/Users/fowler/packages/gemucator/bin/gemucator-run.py", line 14, in <module>
    (locations,bases)=tb_reference_genome.locate_mutation(options.mutation)
  File "/Users/fowler/packages/gemucator/gemucator/core.py", line 127, in locate_mutation
    assert before==bases.translate(), "wildtype amino acid specified in mutation does not match the "+self.genbank_file+" genbank file"
AssertionError: wildtype amino acid specified in mutation does not match the config/H37rV.gbk genbank file
```

Now we can go the other way as well.

```
> gemucator-run.py --location 761153
rpoB_S450
> gemucator-run.py --location 761154
rpoB_S450
> gemucator-run.py --location 761155
rpoB_S450
```

It also handles promoter (nucleotide) mutations. e.g.

```
> gemucator-run.py --mutation pncA_t-12c
pncA_t-12c:
2289252 t
```

Now a single row is returned. Again the code will check that what you give it matches the genbank file! Likewise, it checks that you are giving it a nucleotide.

```
> gemucator-run.py --mutation pncA_x-12c
Traceback (most recent call last):
  File "/Users/fowler/Library/Python/3.5/bin/gemucator-run.py", line 6, in <module>
    exec(compile(open(__file__).read(), __file__, 'exec'))
  File "/Users/fowler/packages/gemucator/bin/gemucator-run.py", line 14, in <module>
    (locations,bases)=tb_reference_genome.locate_mutation(options.mutation)
  File "/Users/fowler/packages/gemucator/gemucator/core.py", line 61, in locate_mutation
    assert before in ['c','t','g','a'], before+" is not a nucleotide!"
AssertionError: x is not a nucleotide!
```
Note that mutation->location->mutation is not uniquely defined for some 'promoters' i.e. the promoter for gene X may lie within the coding region of gene Y which makes assigning it as a CDS mutation or a PROM mutation difficult.

Finally, it will parse insertions and deletions as long as they conform to the format like in the example below.

```
> gemucator-run.py --mutation rpoB_1300_ins_*
rpoB_1300_ins_*:
761106 t
```

This means an insertion (`ins`) of any length (`*`) at nucleotide `1300` in the coding sequence of the `rpoB` gene. You can replace the wildcard with a positive integer to be specific about the number of bases inserted (e.g. for a frame shift). Likewise, for a deletion replace `ins` with `del`.

## Installation

First clone the repository to your local machine

```
> git clone https://github.com/philipwfowler/gemucator.git
```
Now enter the directory and install

```
> cd gemucator
> python3 setup.py install --user
```

The `--user` flag will install the python package in the `$HOME` directory of this user and means you don't need the root password etc. The only dependency is BioPython version 1.70 or newer and the above process will download and install it if it cannot find BioPython on your machine. Now the `gemucator-run.py` script should be in your `$PATH` so try typing one of the examples above!
