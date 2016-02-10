"""
Unit tests for reads objects. This is used for all tests
that can be performed in isolation from input data.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.datamodel.reads as reads
import ga4gh.protocol as protocol


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


class TestSamCigar(unittest.TestCase):
    """
    Test Sam Cigar class handles Cigar mappings correctly

    The integer codes are defined in the SAM spec. Thus, the ordering of
    SamCigar.cigarStrings implicitly implements this spec.
    """

    def testAlignmentMatch(self):
        self.assertEqual(0, reads.SamCigar.ga2int(
            protocol.CigarOperation.ALIGNMENT_MATCH))

        self.assertEqual(protocol.CigarOperation.ALIGNMENT_MATCH,
                         reads.SamCigar.int2ga(0))

    def testInsertion(self):
        self.assertEqual(1, reads.SamCigar.ga2int(
            protocol.CigarOperation.INSERT))

        self.assertEqual(protocol.CigarOperation.INSERT,
                         reads.SamCigar.int2ga(1))

    def testDeletion(self):
        self.assertEqual(2, reads.SamCigar.ga2int(
            protocol.CigarOperation.DELETE))

        self.assertEqual(protocol.CigarOperation.DELETE,
                         reads.SamCigar.int2ga(2))

    def testSkipped(self):
        self.assertEqual(3, reads.SamCigar.ga2int(
            protocol.CigarOperation.SKIP))

        self.assertEqual(protocol.CigarOperation.SKIP,
                         reads.SamCigar.int2ga(3))

    def testSoftClipping(self):
        self.assertEqual(4, reads.SamCigar.ga2int(
            protocol.CigarOperation.CLIP_SOFT))

        self.assertEqual(protocol.CigarOperation.CLIP_SOFT,
                         reads.SamCigar.int2ga(4))

    def testHardClipping(self):
        self.assertEqual(5, reads.SamCigar.ga2int(
            protocol.CigarOperation.CLIP_HARD))

        self.assertEqual(protocol.CigarOperation.CLIP_HARD,
                         reads.SamCigar.int2ga(5))

    def testPadding(self):
        self.assertEqual(6, reads.SamCigar.ga2int(
            protocol.CigarOperation.PAD))

        self.assertEqual(protocol.CigarOperation.PAD,
                         reads.SamCigar.int2ga(6))

    def testSequenceMatch(self):
        self.assertEqual(7, reads.SamCigar.ga2int(
            protocol.CigarOperation.SEQUENCE_MATCH))

        self.assertEqual(protocol.CigarOperation.SEQUENCE_MATCH,
                         reads.SamCigar.int2ga(7))

    def testSequenceMismatch(self):
        self.assertEqual(8, reads.SamCigar.ga2int(
            protocol.CigarOperation.SEQUENCE_MISMATCH))

        self.assertEqual(protocol.CigarOperation.SEQUENCE_MISMATCH,
                         reads.SamCigar.int2ga(8))
