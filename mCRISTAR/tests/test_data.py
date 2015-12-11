## stub for testing example clusters

from os.path import dirname,realpath
from Bio import SeqIO

from ..data import SP04, SP05, SP06, SP10, SP12, SP13, SP14, SP15, SP16, SP17, SP18, SP19, SP21, SP23
from ..data import RBS04, RBS05, RBS06, RBS07, RBS08, RBS09,RBS10,RBS11,RBS12,RBS13,RBS14,RBS15,RBS16,RBS17,RBS18
from ..data import insSP04,insSP05,insSP06,insSP10,insSP12,insSP13,insSP14,insSP15,insSP16,insSP17,insSP18,insSP19,insSP21,insSP23
from ..data import SNP11,SNP12,SNP13,SNP14,SNP15,SNP16,SNP17
from ..data import selective_promoters, nonselective_promoters
from ..data import leu, met
from ..data import crispr_DR, gateway_5prime,gateway_3prime


promoters = [SP04, SP05, SP06, SP10, SP12, SP13, SP14, SP15, SP16, SP17, SP18, SP19, SP21, SP23]
RMS = [RBS04, RBS05, RBS06, RBS07, RBS08, RBS09,RBS10,RBS11,RBS12,RBS13,RBS14,RBS15,RBS16,RBS17,RBS18]
insulators = [insSP04,insSP05,insSP06,insSP10,insSP12,insSP13,insSP14,insSP15,insSP16,insSP17,insSP18,insSP19,insSP21,insSP23]
markers = [leu, met]

def test_selective_sequences():
    """
    test that GB files in data directory are faithful to the manaual sequence in crisprfactor.data.py
    """
    crisprfactordir = dirname(realpath(__file__))
    leucine_gb = crisprfactordir + "/../../data/resistancecassettes/1_LEU2.gb"
    methionine_gb = crisprfactordir + "/../../data/resistancecassettes/2_MET15.gb"
    leu2 = SeqIO.read(leucine_gb,"gb")
    met15 = SeqIO.read(methionine_gb,"gb")
    assert(str(leu2.seq) == str(leu.seq))
    assert(str(met15.seq) == str(met.seq))


def test_sequence_feature_lengths():
    """indivdual features should only have a signel feature,
       selective promoters should have five features,
       nonselective promoters should have six features
    """
    for sequence in promoters + RMS + insulators + markers:
        assert(len(sequence.features) == 1) #single feature
        feature = sequence.features[0]
        featuresize = feature.location.end - feature.location.start
        assert(len(sequence.seq) == featuresize)

        #test the feature-extracted sequence == the whole sequence
        featseq = str(feature.extract(sequence).seq)
        assert(featseq == str(sequence.seq))

    for promoter in selective_promoters:
        assert(len(promoter.features) == 5)

    for promoter in nonselective_promoters:
        assert(len(promoter.features) == 6)


def test_crispr_components():
    #assert full length feature in crispr_DR sequence
    assert(len(crispr_DR.features) == 1)
    assert(len(crispr_DR.features[0].extract(crispr_DR).seq) ==
           len(crispr_DR.seq))

    assert(len(gateway_5prime.features) == 2)
    assert(str(gateway_5prime.features[0].extract(gateway_5prime).seq) == "GGTCTC")
    assert(str(gateway_5prime.features[1].extract(gateway_5prime).seq) == "CCAAAAC")

    assert(len(gateway_3prime.features) == 2)
    assert(str(gateway_3prime.features[0].extract(gateway_3prime).seq) == "GTTTTAGAG")
    assert(str(gateway_3prime.features[1].extract(gateway_3prime).seq) == "GAGACC")



def test_biopython_reverse_complement_adding():
    # test that when adding together sequences with reverse complement you get the correct
    # sequence compositions
    s1 = SP04 + SP04
    s2 = SP04 + SP04.reverse_complement()
    s3 = SP04.reverse_complement() + SP04.reverse_complement()
    assert(str(s1.seq) != str(s2.seq))
    assert(str(s1.seq) == str(s3.reverse_complement().seq))



def test_compound_promoters():
    """
    Promoters consists of a set of features depending on whether or not it is
    selective or non-selective. make sure that when you add them, all sequence
    is accounted for

    """

    #SNP11 = RBS04.reverse_complement() + SP15.reverse_complement() + leu + SP16 + RBS05
    #SNP12 = RBS07.reverse_complement() + SP14.reverse_complement() + insSP14.reverse_complement() + insSP13 + SP13 + RBS08
    #SNP13 = RBS09.reverse_complement() + SP19.reverse_complement() + met + SP12 + RBS10
    #SNP14 = RBS06.reverse_complement() + SP21.reverse_complement() + insSP21.reverse_complement() + insSP23 + SP23 + RBS12
    #SNP15 = RBS13.reverse_complement() + SP17.reverse_complement() + insSP17.reverse_complement() + insSP18 + SP18 + RBS16
    #SNP16 = RBS17.reverse_complement() + SP06.reverse_complement() + insSP06.reverse_complement() + insSP10 + SP10 + RBS18
    #SNP17 = RBS14.reverse_complement() + SP04.reverse_complement() + insSP04.reverse_complement() + insSP05 + SP05 + RBS15

    def check_feature_length(feat, gbk):
        "match feature location values agianst sequences"
        feature_len = feat.location.end - feat.location.start
        featureseq = str(feat.extract(gbk).seq)
        assert(feature_len == len(featureseq))


    for promoter in [SNP11,SNP12,SNP13,SNP14,SNP15,SNP16,SNP17]:
        for feat in promoter.features:
            check_feature_length(feat, gbk=promoter)

    # check SNP11 and leucine
    leu_feat_seq = SNP11.features[2].extract(SNP11)
    assert(str(leu.seq) ==  str(leu_feat_seq.seq))

    # check SNP13 and methionine
    met_feat_seq = SNP13.features[2].extract(SNP13)
    assert(str(met.seq) ==  str(met_feat_seq.seq))

    # #check each feature in a compound promoter
    # RBS07rc, SP14rc ,insSP14rc,insSP13f,SP13f,RBS08f = SNP12.features
    # assert(str(RBS07rc.extract(SNP12).seq) == str(RBS07.reverse_complement().seq))
	#
    # assert(len(insSP14rc.seq) == len(insSP14.seq))
    # assert(str(insSP14rc.seq) == str(insSP14.reverse_complement().seq))








