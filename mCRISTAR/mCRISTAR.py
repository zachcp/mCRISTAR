# -*- coding: utf-8 -*-

#from __future__ import print_function
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from utilities import grouper, fillfeatures, simplegrouper
from data import crisprspacer, selective_promoters, nonselective_promoters, gateway_3prime, gateway_5prime, crispr_DR

class mCRISTAR(object):
    """
    The main CrispFactor object. This object
     1. Designs the CRISP Constructs
            1. Keeps track of the genes
            2. Calculates the intervals to target CRISPR
            3. Checks for Alt Cut Sites on CRISPR
            4. Designs the CRISPR Casettes

     2. Designs the Recombination Cassettes
            1. Find Genes flanking CRISPR Cuts
            2. Design Primer Overhangs for generating the "bridge" priemrs

     3. Outputs a tabular report of CRISPR Sites, Rcombination Casseettes.
        Alt: output Interactive SVG of the process.

    """
    def __init__(self, gbk_file, min_operon_dist, overlaplength=40, maximum_cuts=7, crisprspacer = crisprspacer):
        """

        :return: mCRISTAR Object
        """
        if isinstance(gbk_file, str):
            self.gbk = SeqIO.read(gbk_file,"gb")
        elif isinstance(gbk_file, SeqRecord):
            self.gbk = gbk_file
        else:
            raise(ValueError("gbk value may only be a filename or a SeqRecord"))

        self.overlaplength = overlaplength
        self.min_operon_dist = min_operon_dist
        self.maximum_cuts = maximum_cuts
        self.crisprspacer = crisprspacer
        self.selective_promoters = selective_promoters
        self.nonselective_promoters = nonselective_promoters
        self.gaps = get_gaps(gbk=self.gbk,
                             mindist=self.min_operon_dist,
                             maximum_cuts=self.maximum_cuts)
        self.crisprsites = map(self.get_crispr_sites, self.gaps)
        self.crisprcassetes = self.make_crispr_cassette(crisprsites = self.crisprsites,
                                                        gbk=self.gbk,
                                                        arraysize=4)
        self.bridgeprimers = create_promoter_primers(seqfeatures = self.gaps,
                                                     selective_promoters = self.selective_promoters,
                                                     nonselective_promoters = self.nonselective_promoters,
                                                     gbk = self.gbk,
                                                     overlaplength= self.overlaplength)

        assert(len(self.gaps) == len(self.crisprsites) == len(self.bridgeprimers))
        #assert(len(self.crisprcassetes) ==  (len(self.gaps)/4) + 1)

    def get_crispr_sites(self, gap):
        return find_crispr_site(gap, gbk=self.gbk)

    def make_crispr_cassette(self, crisprsites, gbk, arraysize):
        """
        group CRISPR sequences for use in the cassette. Currently only
        four cuts are possible per cassete with a total of 7 promoter combinations.

        rep-seq1-rep-seq2-rep-seq3-rep-seq4-rep

        """
        crispr_groups = grouper(crisprsites,arraysize)

        crispr_cassettes = []
        for crisprs in crispr_groups:
            # need to obtain the features as follows:
            # <- gateway 5p - {DR - crispr}x - DR - gateway 3p -->
            clist = []
            clist.append(gateway_5prime)

            for idx, crispr in enumerate(crisprs):
                if idx == 0:
                    clist.append(crispr)
                else:
                    clist.append(crispr_DR)
                    clist.append(crispr)

            clist.append(gateway_3prime)

            # combine the list into a single seqrecord
            cassette = reduce(lambda x,y: x+ y, clist)
            # add filler sequences to the crisprcassete
            cassette = fillfeatures(cassette)
            crispr_cassettes.append(cassette)

        return crispr_cassettes


def create_promoter_primers(seqfeatures, selective_promoters, nonselective_promoters, gbk, overlaplength):
    """
    :param seqfeatures: triplet of (prot1,gap,prot2) SequenceFeatures
    :param selective_promoters: Sequence Records of selective promoters
    :param nonselective_promoters: Sequence Records of nonselective promoters
    :param promoterseqs: list of primer adaptors to use the auxotroph genes of interest
    :return:
    """
    def get_selective_promoter():
        for promoter in selective_promoters:
            yield promoter

    def get_nonselective_promoter():
        for promoter in nonselective_promoters:
            yield promoter

    def mod4(x):
        return True if x % 4 == 0 else False

    selective = get_selective_promoter()
    nonselective = get_nonselective_promoter()

    promoter_primers = []
    for i, seqs in enumerate(seqfeatures):
        prot1, gap, prot2 = seqs
        if mod4(i):
            promoter = selective.next()
            primerset = create_primerset(prot1, prot2, promoter, gbk, selection=True)
        else:
            promoter = nonselective.next()
            print promoter
            primerset = create_primerset(prot1, prot2, promoter, gbk, selection=False)

        promoter_primers.append(primerset)

    return promoter_primers


def create_primerset(prot1, prot2, promoter, gbk, overlaplength=40, selection=None):
    """
    Calculate the Primers needed to create the Bridge DNA Using a
    Promoter with a Selective Element.

    :param prot1:
    :param prot2:
    :param promoter:
    :param gbk:
    :param selection:
    :return:
    """

    # promoters either have 6 features in the case of nonselective primers
    # or they have 5 features in the case of selective promoters
    if selection is True:
        rbsF_rc, promoterF_rc, selective, promoterR, rbsR = promoter.features
    else:
        rbsF_rc, promoterF_rc, insF_rc, insR, promoterR, rbsR = promoter.features


    # Get Sequences from the Protein Features USed for Overlaps
    # get final overlaplengths of the 5'->3' prot1
    if prot1.location.strand == 1:
        fprimseq = str(prot1.extract(gbk).seq[-overlaplength:])
    elif prot1.location.strand == -1:
        fprimseq = str(prot1.extract(gbk).reverse_complement().seq[-overlaplength:])
    else:
        raise(ValueError("SeqFeatures Require a Strand"))

    # get first 20 RC'ed bp of prot2
    if prot2.location.strand == 1:
        rprimseq = str(prot2.extract(gbk).reverse_complement().seq[-overlaplength:])
    elif prot2.location.strand == -1:
        rprimseq = str(prot2.extract(gbk).seq[-overlaplength:])
    else:
        raise(ValueError("SeqFeatures Require a Strand"))

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

    if prot1.location.strand == 1 and prot2.location.strand == 1:
        strand_orientation = "sright"
        if selection is True:
            fp = fprimseq +  str(selective.extract(promoter).seq[:18])
            rp = rprimseq +  str(promoter[-18:].seq.reverse_complement())
        else:
            fp = fprimseq +  str(insR.extract(promoter).seq[:18])
            rp = rprimseq +  str(promoter[-18:].seq.reverse_complement())

    elif prot1.location.strand == -1 and prot2.location.strand == -1:
        strand_orientation = "sleft"
        if selection is True:
            fp = fprimseq +  str(promoter[:18].seq)
            rp = rprimseq +  str(selective.extract(promoter)[-18:].seq.reverse_complement())
        else:
            fp = fprimseq +  str(promoter[:18].seq)
            rp = rprimseq +  str(insF_rc.extract(promoter).seq[:18])

    elif prot1.location.strand == -1 and prot2.location.strand == 1:
        strand_orientation = "bi-good"
        fp = fprimseq +  str(promoter[:18].seq)
        rp = rprimseq +  str(promoter[-18:].reverse_complement().seq)

    elif prot1.location.strand == 1 and prot2.location.strand == -1:
        strand_orientation = "bi-bad"
        raise ValueError("gap with inward facing proteins do not need promoters.")
    else:
        raise ValueError("protein must to be oriented")

    return  {"selection": selection,
             "forward": fp,
             "reverse": rp,
             "strandorientation": strand_orientation}

def get_gaps(gbk, mindist):
    """
    return a list of operons as determined by a minimum
    gap distance, mindist, between adjacent proteins

    :param gbk: gbk for finding sites.
    :param mindist: minimum distance between ORFs to look for gaps. Does not apply to beginning/end which require 100 bp.
    :return: list of SeqFeature triplets (prot1, gap, prot2) corresponding to proteins and gaps where
             the gap is bigger than the minimum sequence
    """
    def checkdist(prot1,prot2, mindist):
        "check gap between print1 and prot2 is bigger than minimum"
        return (prot2.location.start - prot1.location.end) > mindist

    def checkorientation(prot1,prot2):
        "ensure orientation is not two proteins pointing toward the gap"
        if (prot1.location.strand == 1) and (prot2.location.strand == -1):
            return False
        else:
            return True

    # keep everything as a seqfeature relative to the top-level info
    cds = [f for f in gbk.features if f.type == "CDS"]
    gaps = []

    for i, _ in enumerate(cds):
        # check for crisprsite upstream of first CDS site
        if i == 0:
            startprotein = cds[i]
            if (startprotein.location.strand == 1) and startprotein.location.start >= 99:
                #create gap triplet by making a dummy protein that is there
                gap = SeqFeature(FeatureLocation(startprotein.location.start-50, startprotein.location.start), type="gap", strand=1)
                dummyprotein = SeqFeature(FeatureLocation(0, startprotein.location.start-50), type="CDS", strand=1)
                gaps.append((dummyprotein, gap, startprotein))

        # check for crisprsite downstream of last CDS site
        if i == len(cds) - 1:
            endprotein = cds[i]
            end_sequence = len(gbk.seq)
            if (endprotein.location.strand == -1) and ((end_sequence - endprotein.location.end) >= 99):
                #create gap triplet by making a dummy protein that is there
                gap = SeqFeature(FeatureLocation(endprotein.location.end, endprotein.location.end + 50), type="gap", strand=1)
                dummyprotein = SeqFeature(FeatureLocation(endprotein.location.end + 50 , endprotein.location.end+100), type="CDS", strand=-1)
                gaps.append((endprotein, gap, dummyprotein))

        # find gaps between two proteins if they are bigger that a minimum site.
        if i > 0:
            prot1 = cds[i-1]
            prot2 = cds[i]
            if checkdist(prot1,prot2,mindist) and checkorientation(prot1,prot2):
                gap = SeqFeature(FeatureLocation(prot1.location.end,prot2.location.start), type="gap", strand=1)
                gaps.append((prot1, gap, prot2))
    return gaps

def find_crispr_site(seqfeats, gbk):
    """

    :param seqfeat: a triplet of SequenceFeatures (prot1, gap, prot2)
    :param gbkfile:
    :return: a SequenceRecord with a dingle feature
    """
    #feature data
    prot1, gap, prot2 = seqfeats
    seqstr   = str(gap.extract(gbk).seq)
    assert(len(seqstr) >= 0)
    seqstart = gap.location.start
    seqend   = gap.location.end

    crisprsites = []
    # forward facing crispr sites
    fpat = re.compile(r'(?=(\S{21}GG))')
    forward_hits = re.finditer(fpat, seqstr)
    for match in forward_hits:
        start, seq = match.start(), match.group(1)

        crisprsite = SeqFeature(FeatureLocation(seqstart + start,
                                                seqstart + start + len(seq) -3 ),
                                type="crisprsite",strand=1)
        crisprsites.append(crisprsite)

    # reverse facing crispr sites
    rpat = re.compile(r'(?=(CC\S{21}))')
    reverse_hits = re.finditer(rpat, seqstr)
    for match in reverse_hits:
        start, seq = match.start(), match.group(1)
        crisprsite = SeqFeature(FeatureLocation(seqstart + start + 3,
                                                seqstart + start + len(seq)),
                                type="crisprsite",strand=-1)
        crisprsites.append(crisprsite)

    #choose a single crispr site closest to the middle
    meanloc = ((seqend - seqstart) / 2 ) + seqstart
    def getdist(x):
        return abs(meanloc - x.location.start)
    dists = map(getdist, crisprsites)
    minindex = dists.index(min(dists))
    crisprsite = crisprsites[minindex]

    #convert crisprsite from a sequence featrue to its own sequence record with a single feature
    crisprsequence = str(crisprsite.extract(gbk).seq)

    return SeqRecord(Seq(crisprsequence),
                     features = [SeqFeature(FeatureLocation(0, len(crisprsequence)), type="crisprsite", strand=1)])





def create_promoter_primers_from_JSON(gaps_from_JSON, selective_promoters, nonselective_promoters, overlaplength):
    """
    :param seqfeatures: json with (prot1,gap,prot2) keys with start end sequence strand keys
    :param selective_promoters: Sequence Records of selective promoters
    :param nonselective_promoters: Sequence Records of nonselective promoters
    :param promoterseqs: list of primer adaptors to use the auxotroph genes of interest
    :return:
    """
    def get_selective_promoter():
        for promoter in selective_promoters:
            yield promoter

    def get_nonselective_promoter():
        for promoter in nonselective_promoters:
            yield promoter

    def mod4(x):
        return True if x % 4 == 0 else False

    selective = get_selective_promoter()
    nonselective = get_nonselective_promoter()

    promoter_primers = []
    for i, gapdata in enumerate(gaps_from_JSON):
        if mod4(i):
            promoter = selective.next()
            primerset = create_primerset_from_JSON(gapdata, promoter, selection=True)
        else:
            promoter = nonselective.next()
            print promoter
            primerset = create_primerset_from_JSON(gapdata, promoter, selection=False)

        promoter_primers.append(primerset)

    return promoter_primers



def create_primerset_from_JSON(gapdata, promoter, overlaplength=40, selection=None):
    """
    Calculate the Primers needed to create the Bridge DNA Using a
    Promoter with a Selective Element.

    :param gapdata: JSON from client
    :param promoter:
    :param selection:
    :return:
    """

    # promoters either have 6 features in the case of nonselective primers
    # or they have 5 features in the case of selective promoters
    if selection is True:
        rbsF_rc, promoterF_rc, selective, promoterR, rbsR = promoter.features
    else:
        rbsF_rc, promoterF_rc, insF_rc, insR, promoterR, rbsR = promoter.features

    #get protein data from the JSON call
    p1 = gapdata['protein1']
    p1_strand = p1['strand']
    p1_sequence = p1['sequence']

    p2 = gapdata['protein2']
    p2_strand = p2['strand']
    p2_sequence = p2['sequence']

    prot1 = SeqRecord(Seq(p1_sequence),
                      features = [SeqFeature(FeatureLocation(0, len(p1_sequence)), type="CDS", strand= p1_strand)])

    prot2 = SeqRecord(Seq(p2_sequence),
                      features = [SeqFeature(FeatureLocation(0, len(p2_sequence)), type="CDS", strand= p2_strand)])


# Get Sequences from the Protein Features USed for Overlaps
    # get final overlaplengths of the 5'->3' prot1
    if prot1.features[0].location.strand == 1:
        fprimseq = str(prot1.seq[-overlaplength:])
    elif prot1.features[0].location.strand == -1:
        fprimseq = str(prot1.reverse_complement().seq[-overlaplength:])
    else:
        raise(ValueError("SeqFeatures Require a Strand"))

    # get first 20 RC'ed bp of prot2
    if prot2.features[0].location.strand == 1:
        rprimseq = str(prot2.reverse_complement().seq[-overlaplength:])
    elif prot2.features[0].location.strand == -1:
        rprimseq = str(prot2.seq[-overlaplength:])
    else:
        raise(ValueError("SeqFeatures Require a Strand"))

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

    if prot1.features[0].location.strand == 1 and prot2.features[0].location.strand == 1:
        strand_orientation = "sright"
        if selection is True:
            fp = fprimseq +  str(selective.extract(promoter).seq[:18])
            rp = rprimseq +  str(promoter[-18:].seq.reverse_complement())
        else:
            fp = fprimseq +  str(insR.extract(promoter).seq[:18])
            rp = rprimseq +  str(promoter[-18:].seq.reverse_complement())

    elif prot1.features[0].location.strand == -1 and prot2.features[0].location.strand == -1:
        strand_orientation = "sleft"
        if selection is True:
            fp = fprimseq +  str(promoter[:18].seq)
            rp = rprimseq +  str(selective.extract(promoter)[-18:].seq.reverse_complement())
        else:
            fp = fprimseq +  str(promoter[:18].seq)
            rp = rprimseq +  str(insF_rc.extract(promoter).seq[:18])

    elif prot1.features[0].location.strand == -1 and prot2.features[0].location.strand == 1:
        strand_orientation = "bi-good"
        fp = fprimseq +  str(promoter[:18].seq)
        rp = rprimseq +  str(promoter[-18:].reverse_complement().seq)

    elif prot1.features[0].location.strand == 1 and prot2.features[0].location.strand == -1:
        strand_orientation = "bi-bad"
        raise ValueError("gap with inward facing proteins do not need promoters.")
    else:
        raise ValueError("protein must to be oriented")

    return  {"promoterid": promoter.id,
            "selection": selection,
             "forward": fp,
             "reverse": rp,
             "strandorientation": strand_orientation}

def find_crispr_site_JSON(JSON_triplet):
    """

    :param seqfeat: a triplet of SequenceFeatures (prot1, gap, prot2)
    :param gbkfile:
    :return: a SequenceRecord with a dingle feature
    """
    #feature data
    seqstr   = JSON_triplet['gap']['sequence']
    seqstart = JSON_triplet['gap']['start']
    seqend   = JSON_triplet['gap']['end']
    assert(len(seqstr) >= 0)


    # keep track of crispr sites with a simple list of dictionaries
    crisprsites = []
    # forward facing crispr sites
    fpat = re.compile(r'(?=(\S{21}GG))')
    forward_hits = re.finditer(fpat, seqstr)

    for match in forward_hits:
        start, seq = match.start(), match.group(1)
        crisprsites.append({"start" : start, "strand": 1, "sequence": seq})

    # reverse facing crispr sites
    rpat = re.compile(r'(?=(CC\S{21}))')
    reverse_hits = re.finditer(rpat, seqstr)
    for match in reverse_hits:
        start, seq = match.start(), match.group(1)
        crisprsites.append({"start" : start, "strand": -1, "sequence": seq})

    #choose a single crispr site closest to the middle
    meanloc = len(seqstr)/2
    def getdist(x):
        return abs(meanloc - x['start'])
    dists = map(getdist, crisprsites)
    minindex = dists.index(min(dists))
    crisprsite = crisprsites[minindex]

    #convert crisprsite from a sequence feature to its own sequence record with a single feature
    if crisprsite['strand'] == 1:
        crisprsequence = crisprsite['sequence'][:-3]
    elif crisprsite['strand'] == -1:
        tempseq = Seq(crisprsite['sequence'][3:])
        crisprsequence = str(tempseq.reverse_complement())
    else:
        raise ValueError("strand must be +1 or -1")

    return SeqRecord(Seq(crisprsequence),
                     features = [SeqFeature(FeatureLocation(0, len(crisprsequence)), type="crisprsite", strand=1)])


def make_crispr_cassette(crisprsites,  arraysize):
    """
    group CRISPR sequences for use in the cassette. Currently only
    four cuts are possible per cassete with a total of 7 promoter combinations.

    rep-seq1-rep-seq2-rep-seq3-rep-seq4-rep

    """
    crispr_groups = simplegrouper(crisprsites)

    crispr_cassettes = []
    for crisprs in crispr_groups:
        # need to obtain the features as follows:
        # <- gateway 5p - {DR - crispr}x - DR - gateway 3p -->
        clist = []
        clist.append(gateway_5prime)

        for idx, crispr in enumerate(crisprs):
            if idx == 0:
                clist.append(crispr)
            else:
                clist.append(crispr_DR)
                clist.append(crispr)

        clist.append(gateway_3prime)

        # combine the list into a single seqrecord
        cassette = reduce(lambda x,y: x+ y, clist)
        # add filler sequences to the crisprcassete
        cassette = fillfeatures(cassette)
        crispr_cassettes.append(cassette)

    return crispr_cassettes


def processGBK(gbk):
    " obtian start/stop/locations for CDS sequences"
    CDSs = [f for f in gbk.features if f.type == "CDS"]

    genes = []
    for idx, cds in enumerate(CDSs):
        start = str(cds.location.start)
        if start == "<0":
            start = 0
        else:
            start = int(start)
        end   = cds.location.end
        strand = cds.location.strand
        genes.append({"id": idx,
                      "start": start,
                      "end": end,
                      "strand": strand,
                      "selected": "true"})
    return genes

def processGaps(gaps, gbk):
    " obtain start/stop/locations for CDS sequences"
    gapdata = []
    for idx, gap in enumerate(gaps):
        prot1, gap, prot2 = gap
        gapdata.append(
        {"id": idx,
         "selected": "true",
         "protein1": {"start": prot1.location.start,
                      "end": prot1.location.end,
                      "sequence": str(prot1.extract(gbk).seq),
                      "strand": prot1.location.strand},
         "gap":     {"start": gap.location.start,
                      "end": gap.location.end,
                      "sequence": str(gap.extract(gbk).seq),
                      "strand": gap.location.strand},
        "protein2": {"start": prot2.location.start,
                      "end": prot2.location.end,
                      "sequence": str(prot2.extract(gbk).seq),
                      "strand": prot2.location.strand}})
    return gapdata


def processCrisprCassettes(cass):
    cassette_data = []
    for feature in cass.features:
        cassette_data.append(
            {"start": 0,
             "end": 0,
             "sequence": str(feature.extract(cass).seq),
             "type": feature.type })
    return cassette_data