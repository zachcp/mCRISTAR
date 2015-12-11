# -*- coding: utf-8 -*-

#from __future__ import print_function
import copy
import re
from itertools import product
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from data import selective_promoters, nonselective_promoters, gateway_3prime, gateway_5prime, crispr_DR


class GBKProcessor(object):
    """
    Handle the processing of the GBK file.
    """
    def __init__(self, gbk):
        """

        :return: mCRISTAR Object
        """
        self.gaps = get_gaps(gbk=gbk, mindist=50)
        self.processedGaps = processGaps(gaps=self.gaps, gbk=gbk)
        self.genes = processGBK(gbk)

    def export(self):
        return {"genes": self.genes,
                "gaps": self.processedGaps}


class GapProcessor(object):
    """
    Obtain a set of gap triplets (prot1, gap, prot2).

    1. Find all of the CRISPR sites within
    each gap. Perform a series of checks to make sure that, when creating the cassettes,
    there are not palindromic sequences added and there are not repeats.

    2. Create the bridge primers for the gaps

    3. Return a json-compatible dicitonary of the cassettes and primers

    """
    def __init__(self, gaptriplets):
        """

        :return: mCRISTAR Object
        """
        self.gaps = gaptriplets
        self.allcrisprsites = [find_crispr_sites_JSON(gap) for gap in self.gaps]
        self.crisprcassettes = self.make_crispr_cassette(self.allcrisprsites)
        self.primers =  create_promoter_primers_from_JSON(self.gaps,
                                                          selective_promoters,
                                                          nonselective_promoters,
                                                          overlaplength=40)

    def make_crispr_cassette(self, crisprsites):
        """
        :param crisprsites: a list of lists of crispr cutsite SeqRecords
        :return: a crispr cassette that has been validated
        """

        #return a group of groups
        cassette_groups = simplegrouper(crisprsites)
        cassette_groups_checked = (choose_crisprsites(cgroup) for cgroup in cassette_groups)
        crispr_cassettes = [make_cassete(crisprsites) for crisprsites in cassette_groups_checked]
        return crispr_cassettes

    def export(self):
        return {"cassettes": [processCrisprCassettes(cass) for cass in self.crisprcassettes],
                "primers": self.primers}



def choose_crisprsites(crisprgroup):
    """
    :param crisprgroup: a list of lists of possible crisprsites that are intended to be joined on a cassette
    :return: tuple of CRISPR sites
    """
    #generator of all possible crispr combinations as a list
    allcrisprcombinations = product(*crisprgroup)

    for candidate in allcrisprcombinations:
        pass_palindromic = confirm_nonpalindromic(candidate)
        pass_extension = confirm_extension(candidate)
        if pass_palindromic and pass_extension:
            return candidate

    raise ValueError("No Acceptable CRISPR combinations were found")

def make_cassete(crisprsites):
    "take validated sites and group them together"
    # need to obtain the features as follows:
    # <- gateway 5p - {DR - crispr}x - DR - gateway 3p -->
    clist = []
    clist.append(gateway_5prime)

    for idx, crispr in enumerate(crisprsites):
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
    return cassette

def confirm_extension(crisprsites):
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


def confirm_nonpalindromic(crisprsites):
    """
    if the end of a crispr site is the RC of the next site,
    the Direct Repeat will be extended. this has been a problem for
    synthesizeing the repeats
    """
    accept = True

    rcdict = {"A":"T",
              "T":"A",
              "C":"G",
              "G":"C"}

    for i,site in enumerate(crisprsites):
        if i > 0:
            seq1 = crisprsites[i-1]
            seq2 = crisprsites[i]
            seq1end = str(seq1[-1])
            seq2start = str(seq2[1])
            if seq1end == rcdict[seq2start]:
                accept = False

    return accept

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

    # Get Sequences from the Protein Features Used for Overlaps
    # get final overlaplengths of the 5'->3' prot1
    # get first 20 RC'ed bp of prot2
    fprimseq = str(prot1.seq[-overlaplength:])
    rprimseq = str(prot2.seq[:overlaplength].reverse_complement())


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

def find_crispr_sites_JSON(JSON_triplet):
    """
    Take a JSON triplet and find all CRISPR sites within each gap

    :param seqfeat: a triplet of SequenceFeatures (prot1, gap, prot2)
    :param gbkfile:
    :return: a list of SequenceRecord with a single SeqFeature
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

    def noBSAIsite(seq):
        if "GGTCTC" in seq:
            return False
        if "GAGACC" in seq:
            return False
        else:
            return True

    def create_crisprsite(crisprsite):
        #convert crisprsite from a sequence feature to its own sequence record with a single feature
        if crisprsite['strand'] == 1:
            crisprsequence = crisprsite['sequence'][:-3]
        elif crisprsite['strand'] == -1:
            tempseq = Seq(crisprsite['sequence'][3:])
            crisprsequence = str(tempseq.reverse_complement())
        else:
            raise ValueError("strand must be +1 or -1")

        # check BSAI site
        if noBSAIsite:
            return SeqRecord(Seq(crisprsequence),
                             features = [SeqFeature(FeatureLocation(0, len(crisprsequence)), type="crisprsite", strand=1)])

    return [create_crisprsite(csite) for csite in crisprsites]


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


def fillfeatures(record, featuretype="filler"):
    """
    Provide a filler feature so that the entire sequence record
    is covered in features. This is useful when you want to iterate over the features a
    and obtain the sequences for each feature without missing a component

    :param gbk: a Biopython SeqRecord
    :param featuretype:  what to set features to
    :return:
    """
    assert(isinstance(record, SeqRecord))

    #record features in local variable
    # remove full length features
    features1 = [feat for feat in record.features if
                (feat.location.start != 0) or (feat.location.end != len(record.seq)) ]

    features = copy.deepcopy(features1)

    # get all features that are not full-length
    feature_coverage = [(feat.location.start, feat.location.end) for feat in features]
    record_length = len(record.seq)

    # calculate new feature intervals
    new_features = []
    number_of_features = len(feature_coverage)

    for idx, (s,e) in enumerate(feature_coverage):
        if idx == 0:
            if s != 0:
                # Add an interval from zero to the start of the first feature if the first feature doesn't start at 0
                new_features.append((0, s))
        elif idx == (number_of_features - 1):
            #add an interval between adjacent features
            new_interval = (feature_coverage[idx-1][1],s)
            new_features.append(new_interval)

            if e != record_length:
                # add an interval from the last feature to the record's end
                new_features.append((e, record_length))

        else:
            new_interval = (feature_coverage[idx-1][1],s)
            new_features.append(new_interval)

    #create new features
    for (fstart,fend) in new_features:
        features.append(
            SeqFeature(FeatureLocation(fstart, fend), type = featuretype, strand=1))

    # sort the features
    features = sorted(features, key=lambda x: x.location.start)

    # assert that all feature components have the same length as the full record
    # i.e we can't handle overlapping features.
    seqs = [str(feat.extract(record).seq) for feat in  features]

    record.features = features
    return record
    # else:
    #     raise ValueError("filling function can only work on GBK files with non-overlapping features.")


def simplegrouper(crisprsites):
    """
    :param crisprsites:
    :param n:
    :return:
    """
    if len(crisprsites) <= 4:
        return [crisprsites]
    elif len(crisprsites) in [5,6]:
        return [crisprsites[0:3],crisprsites[3:]]
    elif len(crisprsites) == 7:
        return [crisprsites[0:4], crisprsites[4:]]
    else:
        raise ValueError("Currently it is only possible to process 7 or fewer crisprsites")


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"


    lists = []
    templist = []
    for idx, it in enumerate(iterable):
        if (idx + 1) % n == 0:
            templist.append(it)
            lists.append(templist)
            templist = []
        else:
            templist.append(it)


    # be sure to include partial lists
    lists.append(templist)

    return [l for l in lists if len(l) > 0 ]