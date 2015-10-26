# from future.standard_library import install_aliases
# install_aliases()
#
# from itertools import zip_longest
import copy
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation



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




# def grouper(iterable, n, fillvalue=None):
#     "Collect data into fixed-length chunks or blocks"
#     # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
#
#     args = [iter(iterable)] * n
#     return zip_longest(*args, fillvalue=fillvalue)

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
