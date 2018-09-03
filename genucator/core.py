#! /usr/bin/env python

import pkg_resources

from Bio import SeqIO

class genucator(object):

    '''The genucator class (short for GenBank Mutation Locator) is designed to tell you the location of a mutation
    in the given GenBank reference genome.

    Args:
        genbank_file (str): The default (and provided) GenBank file is the standard H37rV reference for M. tuberculosis, but any GenBank file should be able to be parsed

    Notes:
        * Mutations are specified in the form gene_mutation e.g. rpoB_S450L. Three forms are accepted: SNP, PROMOTER or INDEL.
        * Amino acids are specified throughout in UPPERCASE and nucleotides in lower case.
        * A SNP mutation specifies the reference and target amino acids and the residue position in the coding sequence e.g. katG_S315T.
        * A PROMOTER mutation requires the reference and target nucleotides to be specified, as well as the location relative to the start codon e.g. inhA_c-15t
        * An INDEL mutation is less specific and requires the location, whether an insertion (ins) or deletion (del) is occuring, and the number of bases (or * as a wildcard for 'any number') e.g. ahpC_66_del_1
        * The code is extremely defensive and if the mutation nomenclature includes the reference amino acid or base, this will be checked against the provided genbank file.
    '''

    def __init__(self, genbank_file=pkg_resources.resource_filename("genucator", "../config/H37rV.gbk")):
        ''' Instantiate an instance of the class by loading the specified GenBank file
        '''

        # use BioPython to load the genbank file
        self.genome=SeqIO.read(genbank_file,'genbank')

        # remember the name of the genbank file for use in assert statements later
        self.genbank_file=genbank_file



    def locate_mutation(self,mutation,nucleotide_mutation=False):

        '''Given a specified mutation, return the location(s) and nucleotide(s) in the reference genome.

        Args:
            mutation (str): a SNP, PROMOTER or INDEL mutation in the format described in the class docstring.
            nucleotide_mutation (bool): if True, forces the code to consider this a simple nucleotide mutation (i.e. not an amino acid). Useful for RNA genes like rrs.
        '''

        # first, parse the mutation
        cols=mutation.split("_")

        # check the mutation at least comprises the expected number of sections
        assert len(cols) in [2,4], "mutation "+mutation+" not in correct format!"

        # the gene/locus name should always be the first component
        gene_name=cols[0]

        # determine if this is a CDS or PROMOTER SNP mutation
        if len(cols)==2:

            before=cols[1][0]
            after=cols[1][-1]
            position=int(cols[1][1:-1])

            # if the position is a negative integer then this is a promoter nucleotide mutation
            if position<0:

                assert before in ['c','t','g','a'], before+" is not a nucleotide!"
                assert after in ['c','t','g','a','*'], after+" is not a nucleotide!"

                mutation_type="PROMOTER"

            # ..otherwise it is an amino acid SNP
            else:

                if nucleotide_mutation:

                    assert before in ['c','t','g','a'], before+" is not a nucleotide!"
                    assert after in ['c','t','g','a','*'], after+" is not a nucleotide!"

                else:

                    assert before in ["!",'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z'], after+" is not an amino acid!"
                    assert after in ['*',"!",'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z'], after+" is not an amino acid!"

                mutation_type="SNP"

        # otherwise it must be an INDEL, which is always nucleotide based
        else:

            position=int(cols[1])

            # be defensive here also!
            assert cols[2] in ["ins","del"], "INDEL must be on the format rpoB_1300_ins_1 i.e. the third element must be ins or del, not "+cols[2]

            if cols[3]!="*":
                assert int(cols[3])>0, "last element in INDEL must be * or a positive integer"

            # only then allow this to be an INDEL!
            mutation_type="INDEL"

        # now that we know what we are looking for, iterate through all the features in the genomes
        for record in self.genome.features:

            # check that the record is a Coding Sequence
            if record.type in ['CDS','rRNA']:

                found_record=False

                # if it is a gene
                if 'gene' in record.qualifiers.keys():
                    if gene_name in record.qualifiers['gene']:
                        found_record=True

                elif 'locus_tag' in record.qualifiers.keys():
                    if gene_name in record.qualifiers['locus_tag']:
                        found_record=True

                if found_record:

                    # the start and end positions in the reference genome
                    start=record.location.start.position
                    end=record.location.end.position

                    # and finally whether the strand (which if -1 means reverse complement)
                    strand=record.location.strand.real

                    # retrieve the coding nucleotides for this gene
                    coding_nucleotides=self.genome[start:end]

                    base_positions=[]

                    if mutation_type=="SNP" and not nucleotide_mutation:

                        if strand==1:

                            base_position=start-1+((3*position)-2)

                            bases=self.genome[base_position:base_position+3].seq

                            if before=="!":
                                assert "*"==bases.translate(), "wildtype amino acid specified in mutation ("+before+") does not match the "+self.genbank_file+" genbank file ("+bases.translate()+")"
                            else:
                                assert before==bases.translate(), "wildtype amino acid specified in mutation ("+before+") does not match the "+self.genbank_file+" genbank file ("+bases.translate()+")"

                            bases=str(bases).lower()

                            for i in bases:
                                base_positions.append(base_position)
                                base_position+=1

                        else:

                            base_position=end+1-((3*position)-2)

                            bases=self.genome[base_position-3:base_position].reverse_complement().seq

                            if before=="!":
                                assert "*"==bases.translate(), "wildtype amino acid specified in mutation ("+before+") does not match the "+self.genbank_file+" genbank file ("+bases.translate()+")"
                            else:
                                assert before==bases.translate(), "wildtype amino acid specified in mutation ("+before+") does not match the "+self.genbank_file+" genbank file ("+bases.translate()+")"

                            bases=str(bases).lower()

                            for i in bases:
                                base_positions.append(base_position)
                                base_position+=-1

                    elif (mutation_type in ["PROMOTER","INDEL"]) or nucleotide_mutation:

                        if strand==1:

                            base_position=start+position

                            bases=self.genome[base_position:base_position+1].seq.lower()

                        else:

                            base_position=end-1-position

                            bases=self.genome[base_position:base_position+1].reverse_complement().seq.lower()

                        # as the reference base is specified for a promoter mutation, be defensive and check it agrees with the genbank file
                        if mutation_type=="PROMOTER":
                            assert bases==before, "wildtype base specified in mutation ("+before+") does not match the "+self.genbank_file+" genbank file ("+bases+")"

                        base_positions.append(base_position)

                    else:

                        raise Exception("mutation "+mutation+" didn't trigger any logic!")

        return(base_positions,bases)

    def identify_gene(self,location):

        assert int(location)>0, "genomic position has to be a positive integer"

        # now that we know what we are looking for, iterate through all the features in the genomes
        for record in self.genome.features:

            # check that the record is a Coding Sequence
            if record.type in ['CDS','rRNA']:

                # the start and end positions in the reference genome
                start=record.location.start.position
                end=record.location.end.position

                if start < location < end:

                    # retrieve the coding nucleotides for this gene
                    coding_nucleotides=self.genome[start:end]

                    gene_name=record.qualifiers['gene'][0]

                    strand=record.location.strand.real

                    if strand==1:

                        position=(location-start+3)//3
                        residue=coding_nucleotides.seq.translate()[position-1]
                    else:
                        position=(location-end-3)//3
                        residue=coding_nucleotides.reverse_complement().seq.translate()[position-1]

                    return(gene_name,residue,position)

        print(location)

        # now that we know what we are looking for, iterate through all the features in the genomes
        for record in self.genome.features:

            # check that the record is a Coding Sequence
            if record.type in ['CDS','rRNA']:

                # the start and end positions in the reference genome
                start=record.location.start.position
                end=record.location.end.position

                strand=record.location.strand.real

                if strand==1:
                    if location>(start-100) and location<start:
                        position=location-start
                        gene_name=record.qualifiers['gene'][0]
                        base=self.genome[location].lower()
                        return(gene_name,base,position)
                else:
                    if location<(end+100) and location>end:
                        position=end-location
                        gene_name=record.qualifiers['gene'][0]
                        base=self.genome[location].lower()
                        return(gene_name,base,position)
