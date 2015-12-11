#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_crisprfactor
----------------------------------

Tests for `crisprfactor` module.
"""

import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from mCRISTAR.data import crisprspacer, selective_promoters, nonselective_promoters, SNP11, SNP12
from mCRISTAR.mCRISTAR import mCRISTAR, get_gaps, find_crispr_site, create_promoter_primers, create_primerset

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
tetaramycin_file = fixname("../data/clusters/Tetaramycin.gb")
def test_testdata():
    "make sure our test data is available"
    gbk = SeqIO.read(open(tetaramycin_file,'r'),"gb")
    assert( isinstance(gbk, SeqRecord))

gbk = SeqIO.read(open(tetaramycin_file,'r'),"gb")
cf = mCRISTAR(gbk, min_operon_dist=150)


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

    cf = mCRISTAR(gbk, min_operon_dist=150)
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
    fakegbk = SeqRecord(Seq("A"*150),
                features = [SeqFeature(FeatureLocation(1, 50), type="CDS", strand=1),
                            SeqFeature(FeatureLocation(100, 150), type="CDS", strand=1)])

    #test selective
    primerset1 = create_primerset(prot1=fakegbk.features[0], prot2=fakegbk.features[1],
                                  gbk=fakegbk, promoter=selective_promoters[0],
                                  selection=True, overlaplength=10)

    expected_primerset1 = {"selection": True,
                           "forward": "A" * 10 + "TCGACGGTCGAGGAGAAC", # beginning of Leu
                           "reverse": "T" * 10 + "TAATTCACCTCCTGAGGC", # end of SNP11
                           "strandorientation": "sright"}

    assert(primerset1.keys() == expected_primerset1.keys())
    assert(primerset1['strandorientation'] == expected_primerset1['strandorientation'])
    assert(primerset1['selection'] == expected_primerset1['selection'])
    assert(primerset1['forward'] == expected_primerset1['forward'])
    assert(primerset1['reverse'] == expected_primerset1['reverse'])

    #test nonselective
    primerset2 = create_primerset(prot1=fakegbk.features[0],prot2=fakegbk.features[1],
                                  gbk=fakegbk,promoter=nonselective_promoters[0],
                                  selection=False,overlaplength=10)

    expected_primerset2 = {"selection": False,
                       "forward": "A" * 10 + "ACCCGGACGCGTGGCACC", #insSP13 sequence from SPN12
                       "reverse": "T" * 10 + "GATCGCCCCTCCTGTGGC", #RC of end of SPN12
                       "strandorientation": "sright"}

    assert(primerset2.keys() == expected_primerset2.keys())
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
    fakegbk = SeqRecord(Seq("A"*150),
                features = [SeqFeature(FeatureLocation(1, 50), type="CDS", strand=-1),
                            SeqFeature(FeatureLocation(100, 150), type="CDS", strand=-1)])

    #test selective
    primerset1 = create_primerset(prot1=fakegbk.features[0],prot2=fakegbk.features[1],
                                  gbk=fakegbk,promoter=selective_promoters[0],
                                  selection=True,overlaplength=10)

    expected_primerset1 = {"selection": True,
                           "forward": "A" * 10 + "TCACAACCCTCCTAGTAA", #begining of SNP11
                           "reverse": "T" * 10 + "TCGACTACGTCGTAAGGC", #RC of the ned of leu2
                           "strandorientation": "sleft"}
    assert(primerset1.keys() == expected_primerset1.keys())
    assert(primerset1['strandorientation'] == expected_primerset1['strandorientation'])
    assert(primerset1['selection'] == expected_primerset1['selection'])
    assert(primerset1['forward'] == expected_primerset1['forward'])
    assert(primerset1['reverse'] == expected_primerset1['reverse'])


    #test nonselective
    primerset2 = create_primerset(prot1=fakegbk.features[0],prot2=fakegbk.features[1],
                                  gbk=fakegbk,promoter=nonselective_promoters[0],
                                  selection=False,overlaplength=10)

    expected_primerset2 = {"selection": False,
                           "forward": "A" * 10 + "AAAGGGTCCTCCTCAACG", #begining of SNP12
                           "reverse": "T" * 10 + "CGCTGATGGCGCTCACCA", #rc of insSP14 sequence from SNP12
                           "strandorientation": "sleft"}

    assert(primerset2.keys() == expected_primerset2.keys())
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
    fakegbk = SeqRecord(Seq("A"*150),
                features = [SeqFeature(FeatureLocation(1, 50), type="CDS", strand=-1),
                            SeqFeature(FeatureLocation(100, 150), type="CDS", strand=1)])

    #test selective
    primerset1 = create_primerset(prot1=fakegbk.features[0],prot2=fakegbk.features[1],
                                  gbk=fakegbk,promoter=selective_promoters[0],
                                  selection=True,overlaplength=10)

    expected_primerset1 = {"selection": True,
                           "forward": "A" * 10 + "TCACAACCCTCCTAGTAA", #begining of SNP11
                           "reverse": "T" * 10 + "TAATTCACCTCCTGAGGC", # end of SNP11
                           "strandorientation": "bi-good"}
    assert(primerset1.keys() == expected_primerset1.keys())
    assert(primerset1['strandorientation'] == expected_primerset1['strandorientation'])
    assert(primerset1['selection'] == expected_primerset1['selection'])
    assert(primerset1['forward'] == expected_primerset1['forward'])
    assert(primerset1['reverse'] == expected_primerset1['reverse'])


    #test nonselective
    primerset2 = create_primerset(prot1=fakegbk.features[0],prot2=fakegbk.features[1],
                                  gbk=fakegbk,promoter=nonselective_promoters[0],
                                  selection=False,overlaplength=10)

    expected_primerset2 = {"selection": False,
                           "forward": "A" * 10 + "AAAGGGTCCTCCTCAACG", #begining of SNP12
                           "reverse": "T" * 10 + "GATCGCCCCTCCTGTGGC", # end of SPN12
                           "strandorientation": "bi-good"}

    assert(primerset2.keys() == expected_primerset2.keys())
    assert(primerset2['strandorientation'] == expected_primerset2['strandorientation'])
    assert(primerset2['selection'] == expected_primerset2['selection'])
    assert(primerset2['forward'] == expected_primerset2['forward'])
    assert(primerset2['reverse'] == expected_primerset2['reverse'])



def test_test_extenion():
    """
    if the beginnings or ends of multiple crispr sites are the same, there will be long
    repeast regions when added next to the Direct Repeat sequence. Check that the first two
    and last two bases of all crispr sites are unique
    """
    site_count = len(crisprsites)
    start_seqs = set([str(site[:2]) for site in crisprsites])
    end_seqs = set([str(site[-2:]) for site in crisprsites])

    if len(start_seqs) != site_count:
        return False
    if len(end_seqs) != site_count:
        return False

    return True

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
