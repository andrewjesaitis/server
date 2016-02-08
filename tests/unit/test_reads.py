"""
Unit tests for reads objects. This is used for all tests
that can be performed in isolation from input data.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.datamodel.reads as reads


class TestParseMalformedBamHeader(unittest.TestCase):
    """
    Unit tests for parsing of malformed bam headers.

    """

    def setUp(self):
        self.headerDictMock = {
            'HD': {'SO': 'coordinate', 'VN': '1.0'},
            'PG': [{'CL': 'bwa index -a bwtsw $reference_fasta',
                    'ID': 'bwa_index',
                    'PN': 'bwa',
                    'VN': '0.5.9-r16'},
                   {'CL': 'bwa aln -q 15 -f $sai_file ' +
                          '$reference_fasta $fastq_file\tPP:bwa_index',
                    'ID': 'bwa_aln_fastq',
                    'PN': 'bwa',
                    'VN': '0.5.9-r16'}],
            'RG': [{'CN': 'WUGSC',
                    'DS': 'SRP001294',
                    'ID': 'SRR062634',
                    'LB': '2845856850',
                    'PI': '206',
                    'PL': 'ILLUMINA',
                    'SM': 'HG00096'}],
            'SQ': [{'LN': 249250621,
                    'M5': '1b22b98cdeb4a9304cb5d48026a85128',
                    'SN': '1',
                    'UR': 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/' +
                          'technical/reference/phase2_reference_assembly' +
                          '_sequence/hs37d5.fa.gz        AS:NCBI37' +
                          '       SP:Human'},
                   {'LN': 243199373,
                    'M5': 'a0d9851da00400dec1098a9255ac712e',
                    'SN': '2',
                    'UR': 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/' +
                          'technical/reference/phase2_reference_assembly' +
                          '_sequence/hs37d5.fa.gz        AS:NCBI37' +
                          '       SP:Human'}]
        }

    def testHeaderLineParsing(self):
        expected = {
            'HD': {'SO': 'coordinate', 'VN': '1.0'},
            'PG': [{'CL': 'bwa index -a bwtsw $reference_fasta',
                    'ID': 'bwa_index',
                    'PN': 'bwa',
                    'VN': '0.5.9-r16'},
                   {'CL': 'bwa aln -q 15 -f $sai_file ' +
                          '$reference_fasta $fastq_file\tPP:bwa_index',
                    'ID': 'bwa_aln_fastq',
                    'PN': 'bwa',
                    'VN': '0.5.9-r16'}],
            'RG': [{'CN': 'WUGSC',
                    'DS': 'SRP001294',
                    'ID': 'SRR062634',
                    'LB': '2845856850',
                    'PI': '206',
                    'PL': 'ILLUMINA',
                    'SM': 'HG00096'}],
            'SQ': [{'LN': 249250621,
                    'M5': '1b22b98cdeb4a9304cb5d48026a85128',
                    'SN': '1',
                    'UR': 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/' +
                          'technical/reference/phase2_reference_assembly' +
                          '_sequence/hs37d5.fa.gz',
                    'AS': 'NCBI37',
                    'SP': 'Human'},
                   {'LN': 243199373,
                    'M5': 'a0d9851da00400dec1098a9255ac712e',
                    'SN': '2',
                    'UR': 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/' +
                          'technical/reference/phase2_reference_assembly' +
                          '_sequence/hs37d5.fa.gz',
                    'AS': 'NCBI37',
                    'SP': 'Human'}]
        }

        for k, v in self.headerDictMock.items():
            if type(v) is dict:
                self.assertEqual(expected[k], reads.parseMalformedBamHeader(v))
            else:
                for idx, item in enumerate(v):
                    self.assertEqual(expected[k][idx],
                                     reads.parseMalformedBamHeader(item))
