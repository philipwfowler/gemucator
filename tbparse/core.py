#! /usr/bin/env python

from Bio import SeqIO


class TBReference(object):
    """docstring for ."""

    def __init__(self, genbank_file="config/H37rV.gbk"):

        self.genome=SeqIO.read(genbank_file,'genbank')


    def genome_position(self,mutation=None,rna_gene=False):

        # first, parse the mutation
        cols=mutation.split("_")

        # check the mutation at least comprises the expected number of sections
        assert len(cols) in [2,4], "mutation "+mutation+" not in correct format!"

        # the gene/locus name should always be the first component
        gene_name=cols[0]

        # if there are only two components it is a SNP (including promoter)
        if len(cols)==2:

            residue_before=cols[1][0]
            residue_after=cols[1][-1]
            residue_position=int(cols[1][1:-1])

            # be defensive!
            assert residue_before in ['*',"!",'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z'], after+" is not an amino acid!"
            assert residue_after in ['*',"!",'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z'], after+" is not an amino acid!"
            assert residue_position>-100, "position must be ≥ -100"

        # iterate through all the features in the genomes
        for record in self.genome.features:

            # check that the record is a Coding Sequence and it is also a gene
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

                    coding_nucleotides=reference_genome[start:end]

                    print(coding_nucleotides.seq)

                    # retrieve the number of the first protein coding codon (usually 1)
                    # if rna_gene:
                    #     first_amino_acid_position=1
                    #     default_promoter_length=0
                    # else:
                    #     codon_start=record.qualifiers['codon_start']
                    #     first_amino_acid_position=int(codon_start[0])
                    #     default_promoter_length=100
