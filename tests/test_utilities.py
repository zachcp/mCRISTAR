"""
This file tests the utilities namespace

"""
from mCRISTAR.utilities import grouper, fillfeatures, simplegrouper
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from test_mCRISTAR import gbk, cf


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
    assert(grouper(l1) == [[1]])
    assert(grouper(l2) == [[1,2]])
    assert(grouper(l3) == [[1,2,3]])
    assert(grouper(l4) == [[1,2,3,4]])
    assert(grouper(l5) == [[1,2,3],[4,5]])
    assert(grouper(l6) == [[1,2,3], [4,5,6]])
    assert(grouper(l7) == [[1,2,3,4],[5,6,7,8]])


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
