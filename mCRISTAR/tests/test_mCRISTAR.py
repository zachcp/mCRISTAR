#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_crisprfactor
----------------------------------

Tests for `crisprfactor` module.
"""

import os
import re
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from ..data import crisprspacer, selective_promoters, nonselective_promoters, SNP11, SNP12
from ..mCRISTAR import GBKProcessor, GapProcessor, get_gaps, grouper, simplegrouper, confirm_extension, fillfeatures, create_primerset_from_JSON

def fixname(filename):
    "alter a filename to correspond to the test directory"
    currentdir = os.path.dirname(os.path.realpath(__file__))
    return currentdir + "/" + filename


################################################################################
################################################################################
# global data
#
#
#
#
#
#
# Load Test Datafile
tetaramycin_file = fixname("../../data/clusters/Tetaramycin.gb")
def test_testdata():
    "make sure our test data is available"
    gbk = SeqIO.read(open(tetaramycin_file,'r'),"gb")
    assert( isinstance(gbk, SeqRecord))

gbk = SeqIO.read(open(tetaramycin_file,'r'),"gb")
gbkProcessor = GBKProcessor(gbk)


def test_GBKProcessor():
    "Basic test to load a GBK file and test the shape of the output"
    data = gbkProcessor.export()
    genedata = data['genes']
    gapdata = data['gaps']

    assert(set(data.keys()) == set(['genes','gaps']))
    assert(set(genedata[0].keys()) == set(['start', 'selected', 'end', 'id', 'strand']))
    assert(set(gapdata[0].keys()) == set(['selected', 'protein2', 'id', 'protein1', 'gap']))


def test_basic_GapProcessor():
    "Test Basic GapPRocessor Input and Data Structure"
    # if you submit greater than 7 gaps a value error it raised
    with pytest.raises(ValueError):
        GapProcessor(gbkProcessor.export()['gaps'])

    gapProcessor = GapProcessor(gbkProcessor.export()['gaps'][:7])
    gpExport = gapProcessor.export()
    assert(set(gpExport.keys()) == set(['primers', 'cassettes']))
    #primers are list of dictionaries
    assert(set(gpExport['primers'][0].keys()) == set(['forward', 'strandorientation', 'selection', 'reverse', 'promoterid']))
    #cassettes are list of lists of dictionaries
    assert(set(gpExport['cassettes'][0][0].keys()) == set(['start', 'type', 'end', 'sequence']))


################################################################################
################################################################################
# Tests
#
#
#
#
#
#
#
#
#  Tests Support Functions
def test_get_gap():
    "make sure gaps are the correct size, return SeqFeatures, and are the correct distance apart"
    gaps_150 = get_gaps(gbk, 150, maximum_cuts=100)
    gaps_200 = get_gaps(gbk, 200, maximum_cuts=100)

    # check length
    assert(len(gaps_150)== 6)
    assert(len(gaps_200)== 4)

    #check return classes
    for (p1, gap, p2) in gaps_150:
        assert(isinstance(p1, SeqFeature))
        assert(isinstance(gap, SeqFeature))
        assert(isinstance(p2, SeqFeature))
        assert(p2.location.start - p1.location.end >= 150)
        assert(gap.location.start == p1.location.end)
        assert(gap.location.end == p2.location.start)

    #check return classes
    for (p1, gap, p2) in gaps_200:
        assert(isinstance(p1, SeqFeature))
        assert(isinstance(gap, SeqFeature))
        assert(isinstance(p2, SeqFeature))
        assert(p2.location.start - p1.location.end >= 200)
        assert(gap.location.start == p1.location.end)
        assert(gap.location.end == p2.location.start)

def test_CRISPR_regex():
    """
    Test the Crispr Finding logic.

    Regex Forwad and Reverse
    """
    fcrisprs = "TAAAAAAAAAAAAAAAAAAAAAAAGGTTTTTTTTTTTTTTTTTTTTTTTTGGCCCCCCCCCCCCCCCCCCCCCCCGG"
    fpat = re.compile(r'(?=(\S{21}GG))')
    fmatches = re.finditer(fpat, fcrisprs)
    fmatches = [m.group(1) for m in fmatches]
    assert(len(fmatches) == 3)
    assert(fmatches[0] == "AAAAAAAAAAAAAAAAAAAAAGG")
    assert(fmatches[1] == "TTTTTTTTTTTTTTTTTTTTTGG")
    assert(fmatches[2] == "CCCCCCCCCCCCCCCCCCCCCGG")

    rcrisprs = "CCCAAAAAAAAAAAAAAAAAAAAAAAACCTTTTTTTTTTTTTTTTTTTTTTTTCCTT"
    rpat = re.compile(r'(?=(CC\S{21}))')
    rmatches = re.finditer(rpat, rcrisprs)
    rmatches = [m.group(1) for m in rmatches]
    assert(len(rmatches) == 3)
    assert(rmatches[0] == "CCCAAAAAAAAAAAAAAAAAAAA")
    assert(rmatches[1] == "CCAAAAAAAAAAAAAAAAAAAAA")
    assert(rmatches[2] == "CCTTTTTTTTTTTTTTTTTTTTT")


def test_crispr_cassettes():
    """
    test the crispcassettes of the final data
    """

    cf = GBKProcessor(gbk, min_operon_dist=150)

    cassettes = cf.crisprcassetes

    for cassette in cassettes:
        features = cassette.features
        feature_seqs = [str(f.extract(cassette).seq) for f in features]
        #assert that features, when added together are the same as the full feature
        assert("".join(feature_seqs) == str(cassette.seq))

        #assert two BSA sites, two cut sites, and a minimum of one crispr cut site
        assert(len([f for f in features if f.type == "BSAI Recognition"]) == 2)
        assert(len([f for f in features if f.type == "BSAI Cut Site"]) == 2)
        assert(len([f for f in features if f.type == "crisprsite"]) >= 1)


def test_promoter_primers():
    pass



def test_create_primerset_colinear_right():
    """

    # Get Sequences from the Insert
    #  ---p1--->  ----p2---> Single, right-facing sright
    #  <--p1----  <---p2---- Single,  left-facing sleft
    #  <--p1----  ----p2---> Bidirectional.       bi-good
    #  ---p1--->  <---p2---- Don't Use.           bi-bad
    #
    # Promoter Structure: Selective
    # rbsF_rc  promoterF_rc selective promoterR  rbsR
    # <------  <----------- --------> -------->  ----->
    #
    # Promoter Structure: Nonselective
    # rbsF_rc  promoterF_rc insF_rc insR  promoterR  rbsR
    # <------  <----------- <------ ----> -------->  ----->


    """
    # fakegbk = SeqRecord(Seq("A"*150),
    #             features = [SeqFeature(FeatureLocation(1, 50), type="CDS", strand=1),
    #                         SeqFeature(FeatureLocation(100, 150), type="CDS", strand=1)])
    fakegapdata = {'protein1': {'start': 1,
                        'end': 50,
                        'strand': 1,
                        'sequence': "A"*50},
           'gap':      {'start': 51,
                        'end': 99,
                        'strand': 1,
                        'sequence': "A"*50},
           'protein2': {'start': 100,
                        'end': 150,
                        'strand': 1,
                        'sequence': "A"*50},
           'id': 0,
           'selected': 'true'}
    #test selective
    # primerset1 = create_primerset(prot1=fakegbk.features[0], prot2=fakegbk.features[1],
    #                               gbk=fakegbk, promoter=selective_promoters[0],
    #                               selection=True, overlaplength=10)
    primerset1 = create_primerset_from_JSON(gapdata = fakegapdata,
                                    promoter = selective_promoters[0],
                                    selection=True,
                                    overlaplength=10)

    expected_primerset1 = {"selection": True,
                           "forward": "A" * 10 + "TCGACGGTCGAGGAGAAC", # beginning of Leu
                           "reverse": "T" * 10 + "TAATTCACCTCCTGAGGC", # end of SNP11
                           "promoterid": "SNP11",
                           "strandorientation": "sright"}

    assert(set(primerset1.keys()) == set(expected_primerset1.keys()))
    assert(primerset1['strandorientation'] == expected_primerset1['strandorientation'])
    assert(primerset1['selection'] == expected_primerset1['selection'])
    assert(primerset1['forward'] == expected_primerset1['forward'])
    assert(primerset1['reverse'] == expected_primerset1['reverse'])

    #test nonselective
    # primerset2 = create_primerset(prot1=fakegbk.features[0],prot2=fakegbk.features[1],
    #                               gbk=fakegbk,promoter=nonselective_promoters[0],
    #                               selection=False,overlaplength=10)

    primerset2 = create_primerset_from_JSON(gapdata = fakegapdata,
                                    promoter = nonselective_promoters[0],
                                    selection=False,
                                    overlaplength=10)

    expected_primerset2 = {"selection": False,
                       "forward": "A" * 10 + "ACCCGGACGCGTGGCACC", #insSP13 sequence from SPN12
                       "reverse": "T" * 10 + "GATCGCCCCTCCTGTGGC", #RC of end of SPN12
                       "promoterid": "SNP11",
                       "strandorientation": "sright"}

    assert(set(primerset2.keys()) == set(expected_primerset2.keys()))
    assert(primerset2['strandorientation'] == expected_primerset2['strandorientation'])
    assert(primerset2['selection'] == expected_primerset2['selection'])
    assert(primerset2['forward'] == expected_primerset2['forward'])
    assert(primerset2['reverse'] == expected_primerset2['reverse'])



def test_create_primerset_colinear_left():
    """
    # Get Sequences from the Insert
    #  ---p1--->  ----p2---> Single, right-facing sright
    #  <--p1----  <---p2---- Single,  left-facing sleft
    #  <--p1----  ----p2---> Bidirectional.       bi-good
    #  ---p1--->  <---p2---- Don't Use.           bi-bad
    #
    # Promoter Structure: Selective
    # rbsF_rc  promoterF_rc selective promoterR  rbsR
    # <------  <----------- --------> -------->  ----->
    #
    # Promoter Structure: Nonselective
    # rbsF_rc  promoterF_rc insF_rc insR  promoterR  rbsR
    # <------  <----------- <------ ----> -------->  ----->


    """
    fakegapdata = {'protein1': {'start': 1,
                            'end': 50,
                            'strand': -1,
                            'sequence': "A"*50},
               'gap':      {'start': 51,
                            'end': 99,
                            'strand': 1,
                            'sequence': "A"*50},
               'protein2': {'start': 100,
                            'end': 150,
                            'strand': -1,
                            'sequence': "A"*50},
               'id': 0,
               'selected': 'true'}
    # fakegbk = SeqRecord(Seq("A"*150),
    #             features = [SeqFeature(FeatureLocation(1, 50), type="CDS", strand=-1),
    #                         SeqFeature(FeatureLocation(100, 150), type="CDS", strand=-1)])

    #test selective
    # primerset1 = create_primerset(prot1=fakegbk.features[0],prot2=fakegbk.features[1],
    #                               gbk=fakegbk,promoter=selective_promoters[0],
    #                               selection=True,overlaplength=10)
    primerset1 = create_primerset_from_JSON(gapdata = fakegapdata,
                                        promoter = selective_promoters[0],
                                        selection=True,
                                        overlaplength=10)

    expected_primerset1 = {"selection": True,
                           "forward": "A" * 10 + "TCACAACCCTCCTAGTAA", #begining of SNP11
                           "reverse": "T" * 10 + "TCGACTACGTCGTAAGGC", #RC of the ned of leu2
                           "promoterid": "SNP11",
                           "strandorientation": "sleft"}

    assert(set(primerset1.keys()) == set(expected_primerset1.keys()))
    assert(primerset1['strandorientation'] == expected_primerset1['strandorientation'])
    assert(primerset1['selection'] == expected_primerset1['selection'])
    assert(primerset1['forward'] == expected_primerset1['forward'])
    assert(primerset1['reverse'] == expected_primerset1['reverse'])


    #test nonselective
    # primerset2 = create_primerset(prot1=fakegbk.features[0],prot2=fakegbk.features[1],
    #                               gbk=fakegbk,promoter=nonselective_promoters[0],
    #                               selection=False,overlaplength=10)

    primerset2 = create_primerset_from_JSON(gapdata = fakegapdata,
                                        promoter = nonselective_promoters[0],
                                        selection=False,
                                        overlaplength=10)

    expected_primerset2 = {"selection": False,
                           "forward": "A" * 10 + "AAAGGGTCCTCCTCAACG", #begining of SNP12
                           "reverse": "T" * 10 + "CGCTGATGGCGCTCACCA", #rc of insSP14 sequence from SNP12
                           "promoterid": "SNP12",
                           "strandorientation": "sleft"}

    assert(set(primerset2.keys()) == set(expected_primerset2.keys()))
    assert(primerset2['strandorientation'] == expected_primerset2['strandorientation'])
    assert(primerset2['selection'] == expected_primerset2['selection'])
    assert(primerset2['forward'] == expected_primerset2['forward'])
    assert(primerset2['reverse'] == expected_primerset2['reverse'])


def test_create_primerset_colinear_bidirectional():
    """
    # Get Sequences from the Insert
    #  ---p1--->  ----p2---> Single, right-facing sright
    #  <--p1----  <---p2---- Single,  left-facing sleft
    #  <--p1----  ----p2---> Bidirectional.       bi-good
    #  ---p1--->  <---p2---- Don't Use.           bi-bad
    #
    # Promoter Structure: Selective
    # rbsF_rc  promoterF_rc selective promoterR  rbsR
    # <------  <----------- --------> -------->  ----->
    #
    # Promoter Structure: Nonselective
    # rbsF_rc  promoterF_rc insF_rc insR  promoterR  rbsR
    # <------  <----------- <------ ----> -------->  ----->

    """
    fakegapdata = {'protein1': {'start': 1,
                                'end': 50,
                                'strand': -1,
                                'sequence': "A"*50},
                   'gap':      {'start': 51,
                                'end': 99,
                                'strand': 1,
                                'sequence': "A"*50},
                   'protein2': {'start': 100,
                                'end': 150,
                                'strand': 1,
                                'sequence': "A"*50},
                   'id': 0,
                   'selected': 'true'}

    primerset1 = create_primerset_from_JSON(gapdata = fakegapdata,
                                            promoter = selective_promoters[0],
                                            selection=True,
                                            overlaplength=10)

    expected_primerset1 = {"selection": True,
                           "forward": "A" * 10 + "TCACAACCCTCCTAGTAA", #begining of SNP11
                           "reverse": "T" * 10 + "TAATTCACCTCCTGAGGC", # end of SNP11
                           "promoterid": "SNP11",
                           "strandorientation": "bi-good"}

    assert(set(primerset1.keys()) == set(expected_primerset1.keys()))
    assert(primerset1['strandorientation'] == expected_primerset1['strandorientation'])
    assert(primerset1['selection'] == expected_primerset1['selection'])
    assert(primerset1['forward'] == expected_primerset1['forward'])
    assert(primerset1['reverse'] == expected_primerset1['reverse'])


    #test nonselective
    primerset2 = create_primerset_from_JSON(gapdata = fakegapdata,
                                            promoter = nonselective_promoters[0],
                                            selection=False,
                                            overlaplength=10)

    expected_primerset2 = {"selection": False,
                           "forward": "A" * 10 + "AAAGGGTCCTCCTCAACG", #begining of SNP12
                           "reverse": "T" * 10 + "GATCGCCCCTCCTGTGGC", # end of SPN12
                           "promoterid": "SNP12",
                           "strandorientation": "bi-good"}

    assert(set(primerset2.keys()) == set(expected_primerset2.keys()))
    assert(primerset2['strandorientation'] == expected_primerset2['strandorientation'])
    assert(primerset2['selection'] == expected_primerset2['selection'])
    assert(primerset2['forward'] == expected_primerset2['forward'])
    assert(primerset2['reverse'] == expected_primerset2['reverse'])



def test_confirm_extension():
    """
    if the beginnings or ends of multiple crispr sites are the same, there will be long
    repeast regions when added next to the Direct Repeat sequence. Check that the first two
    and last two bases of all crispr sites are unique
    """
    testlist1 = ["ACT","ACT","ACT"] # same seqs
    testlist2 = ["ACT","CAT","GGG"] # different sets of first two and last two letters
    testlist3 =  ["ACT","ACG","GTG"] # same first two letters of list
    testlist4 =  ["ACT","GCT","GTG"] # same last two letters of list

    assert(confirm_extension(testlist1) == False)
    assert(confirm_extension(testlist2) == True)
    assert(confirm_extension(testlist3) == False)
    assert(confirm_extension(testlist4) == False)

# def test_find_crispr_site():
#     """
#     be sure that crispr sites are the correct size and they are adjacent to a PAM sequence
#     """
#
#     crisprsites = cf.crisprsites
#
#     #check that the CRISPR sites are seqrecords
#     for c in crisprsites:
#         assert(isinstance(c, SeqRecord))
#         assert(len(c.seq) == 20)
#
#     # find the {primer}NGG or CCN{primer} for each crisprsite
#     for cs in crisprsites:
#         if cs.location.strand == 1:
#             #expect PAM site (NGG) downstream
#             pamstart = cs.location.end + 1
#             pamend   = cs.location.end + 3
#             tempfeature = SeqFeature(FeatureLocation(pamstart,pamend), type="pamsite",strand=1)
#             tempseq = str(tempfeature.extract(cf.gbk).seq)
#             assert(tempseq == "GG")
#         else:
#             # expect PAM site (CCN) upstream of the sequence
#             pamstart = cs.location.start - 3
#             pamend   = cs.location.start - 1
#             tempfeature = SeqFeature(FeatureLocation(pamstart,pamend), type="pamsite",strand=1)
#             tempseq = str(tempfeature.extract(cf.gbk).seq)
#             assert(tempseq == "CC")


def test_simplegrouper():
    """ test that the grouper function groups to ther correct level and
     has the expected behcavior for the final group
    """
    l1 = [1]
    l2 = [1,2]
    l3 = [1,2,3]
    l4 = [1,2,3,4]
    l5 = [1,2,3,4,5]
    l6 = [1,2,3,4,5,6]
    l7 = [1,2,3,4,5,6,7]
    assert(simplegrouper(l1) == [[1]])
    assert(simplegrouper(l2) == [[1,2]])
    assert(simplegrouper(l3) == [[1,2,3]])
    assert(simplegrouper(l4) == [[1,2,3,4]])
    assert(simplegrouper(l5) == [[1,2,3],[4,5]])
    assert(simplegrouper(l6) == [[1,2,3], [4,5,6]])
    assert(simplegrouper(l7) == [[1,2,3,4],[5,6,7]])


def test_simplegrouper_lists():
    """ test that the grouper function groups to ther correct level and
     has the expected behcavior for the final group
    """
    l1 = [[1]]
    l2 = [[1],[2]]
    l3 = [[1],[2],[3]]
    l4 = [[1],[2],[3],[4]]
    l5 = [[1],[2],[3],[4],[5]]
    l6 = [[1],[2],[3],[4],[5],[6]]
    l7 = [[1],[2],[3],[4],[5],[6],[7]]
    assert(simplegrouper(l1) == [[[1]]])
    assert(simplegrouper(l2) == [[[1],[2]]])
    assert(simplegrouper(l3) == [[[1],[2],[3]]])
    assert(simplegrouper(l4) == [[[1],[2],[3],[4]]])
    assert(simplegrouper(l5) == [[[1],[2],[3]],[[4],[5]]])
    assert(simplegrouper(l6) == [[[1],[2],[3]], [[4],[5],[6]]])
    assert(simplegrouper(l7) == [[[1],[2],[3],[4]],[[5],[6],[7]]])



def test_grouper():
    """ test that the grouper function groups to ther correct level and
     has the expected behcavior for the final group
    """
    l1 = [1,2,3,4,5,6,7,8]
    assert(grouper(l1,4) == [[1,2,3,4],[5,6,7,8]])
    assert(grouper(l1,3) == [[1,2,3],[4,5,6],[7,8]])


def test_fillfeatures():
    """
    test that the feature filling behavior is correct.

    """
    testrecord = SeqRecord(Seq("A"*10),
                           features= [SeqFeature(FeatureLocation(0,4), type="CDS", strand=1),
                                      SeqFeature(FeatureLocation(8,10), type="CDS",strand=1)])

    assert(len(testrecord.features) == 2)

    #test basic  addition
    newrecord = fillfeatures(testrecord)
    assert(len(newrecord.features) == 3)

    #extract sequences and assert they are the same as the full length
    seqs = [feat.extract(newrecord) for feat in newrecord.features]
    seqsum = sum(len(s.seq) for s in seqs)
    assert(seqsum == len(newrecord.seq))

    #assert that extracted feature sequences, when added together are the same as the original sequence
    seqs_text = [str(s.seq) for s in seqs]
    assert("".join(seqs_text) == str(newrecord.seq))


# def test_fillfeatures_on_crisprcassettes():
#     """
#     use the tetaramycin gene clusters as a test
#     """
#     for cassette in cf.crisprcassetes:
#
#     newrecord = fillfeatures(gbk)
#
#     #extract sequences and assert they are the same as the full length
#     seqs = [feat.extract(newrecord) for feat in newrecord.features]
#     seqsum = sum(len(s.seq) for s in seqs)
#     assert(seqsum == len(newrecord.seq))
#
#     #assert that extracted feature sequences, when added together are the same as the original sequence
#     seqs_text = [str(s.seq) for s in seqs]
#     assert("".join(seqs_text) == str(newrecord.seq))
