# genucator

genucator is short for "Genbank Mutation Locator". It is a simple Python3 class for incorporating into tools that, if you give it a mutation, it will tell its location in the reference genome. The `genucator` class accepts a path to a genbank file; since I am working with /M. tuberculosis/ this is the H37rV genbank file by default, but any genbank file should work.

The package comes with a simple script called `genucator-run.py` that shows how it works. All these examples are for TB.

First, let's see what happens when we give it an amino acid mutation (which has to be, be definition, in the coding sequence of a gene).

```
> genucator-run.py --mutation rpoB_S450L
rpoB_S450L:
761153 t
761154 c
761155 g
```
It returns three rows, since there are three bases in the triplet, each with the position in the H37rV reference genome and the reference base. Note that the code is very defensive and checks that `tcg` is a Serine, which happily it is. If we get the mutation wrong, then the code will stop and catch fire.

```
> genucator-run.py --mutation rpoB_K450L
Traceback (most recent call last):
  File "/Users/fowler/Library/Python/3.5/bin/genucator-run.py", line 6, in <module>
    exec(compile(open(__file__).read(), __file__, 'exec'))
  File "/Users/fowler/packages/genucator/bin/genucator-run.py", line 14, in <module>
    (locations,bases)=tb_reference_genome.locate_mutation(options.mutation)
  File "/Users/fowler/packages/genucator/genucator/core.py", line 127, in locate_mutation
    assert before==bases.translate(), "wildtype amino acid specified in mutation does not match the "+self.genbank_file+" genbank file"
AssertionError: wildtype amino acid specified in mutation does not match the config/H37rV.gbk genbank file
```

It also handles promoter (nucleotide) mutations. e.g.

```
> genucator-run.py --mutation pncA_t-12c
pncA_t-12c:
2289252 t
```

Now a single row is returned. Again the code will check that what you give it matches the genbank file! Likewise, it checks that you are giving it a nucleotide.

```
> genucator-run.py --mutation pncA_x-12c
Traceback (most recent call last):
  File "/Users/fowler/Library/Python/3.5/bin/genucator-run.py", line 6, in <module>
    exec(compile(open(__file__).read(), __file__, 'exec'))
  File "/Users/fowler/packages/genucator/bin/genucator-run.py", line 14, in <module>
    (locations,bases)=tb_reference_genome.locate_mutation(options.mutation)
  File "/Users/fowler/packages/genucator/genucator/core.py", line 61, in locate_mutation
    assert before in ['c','t','g','a'], before+" is not a nucleotide!"
AssertionError: x is not a nucleotide!
```
Finally, it will parse insertions and deletions as long as they conform to the format like in the example below.

```
> genucator-run.py --mutation rpoB_1300_ins_*
rpoB_1300_ins_*:
761106 t
```

This means an insertion (`ins`) of any length (`*`) at nucleotide `1300` in the coding sequence of the `rpoB` gene. You can replace the wildcard with a positive integer to be specific about the number of bases inserted (e.g. for a frame shift). Likewise, for a deletion replace `ins` with `del`. 

## Installation

First clone the repository to your local machine

```
> git clone https://github.com/philipwfowler/genucator.git
```
Now enter the directory and install 

```
> cd genucator
> python3 setup.py install --user
```

The `--user` flag will install the python package in the `$HOME` directory of this user and means you don't need the root password etc. The only dependency is BioPython version 1.70 or newer and the above process will download and install it if it cannot find BioPython on your machine. Now the `genucator-run.py` script should be in your `$PATH` so try typing one of the examples above!




