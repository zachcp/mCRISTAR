# Data File with information needed for cloning
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

# CRISPR Info
crisprspacer = "GTTTTAGAGCTATGCTGTTTTGAATGGTCCCAAAAC"

################################################################################################
################################################################################################
#
# Selective yeast Markers
#
#
#create a new Sequence record with a single feature and a single sequence
leu = SeqRecord(Seq("TCGACGGTCGAGGAGAACTTCTAGTATATCCACATACCTAATATTATTGCCTTATTAAAAATGGAATCCCAACAATTACATCAAAATCCACATTCTCTTCAAAATCAATTGTCCTGTACTTCCTTGTTCATGTGTGTTCAAAAACGTTATATTTATAGGATAATTATACTCTATTTCTCAACAAGTAATTGGTTGTTTGGCCGAGCGGTCTAAGGCGCCTGATTCAAGAAATATCTTGACCGCAGTTAACTGTGGGAATACTCAGGTATCGTAAGATGCAAGAGTTCGAATCTCTTAGCAACCATTATTTTTTTCCTCAACATAACGAGAACACACAGGGGCGCTATCGCACAGAATCAAATTCGATGACTGGAAATTTTTTGTTAATTTCAGAGGTCGCCTGACGCATATACCTTTTTCAACTGAAAAATTGGGAGAAAAAGGAAAGGTGAGAGGCCGGAACCGGCTTTTCATATAGAATAGAGAAGCGTTCATGACTAAATGCTTGCATCACAATACTTGAAGTTGACAATATTATTTAAGGACCTATTGTTTTTTCCAATAGGTGGTTAGCAATCGTCTTACTTTCTAACTTTTCTTACCTTTTACATTTCAGCAATATATATATATATTTCAAGGATATACCATTCTAATGTCTGCCCCTATGTCTGCCCCTAAGAAGATCGTCGTTTTGCCAGGTGACCACGTTGGTCAAGAAATCACAGCCGAAGCCATTAAGGTTCTTAAAGCTATTTCTGATGTTCGTTCCAATGTCAAGTTCGATTTCGAAAATCATTTAATTGGTGGTGCTGCTATCGATGCTACAGGTGTCCCACTTCCAGATGAGGCGCTGGAAGCCTCCAAGAAGGTTGATGCCGTTTTGTTAGGTGCTGTGGCTGGTCCTAAATGGGGTACCGGTAGTGTTAGACCTGAACAAGGTTTACTAAAAATCCGTAAAGAACTTCAATTGTACGCCAACTTAAGACCATGTAACTTTGCATCCGACTCTCTTTTAGACTTATCTCCAATCAAGCCACAATTTGCTAAAGGTACTGACTTCGTTGTTGTCAGAGAATTAGTGGGAGGTATTTACTTTGGTAAGAGAAAGGAAGACGATGGTGATGGTGTCGCTTGGGATAGTGAACAATACACCGTTCCAGAAGTGCAAAGAATCACAAGAATGGCCGCTTTCATGGCCCTACAACATGAGCCACCATTGCCTATTTGGTCCTTGGATAAAGCTAATCTTTTGGCCTCTTCAAGATTATGGAGAAAAACTGTGGAGGAAACCATCAAGAACGAATTCCCTACATTGAAGGTTCAACATCAATTGATTGATTCTGCCGCCATGATCCTAGTTAAGAACCCAACCCACCTAAATGGTATTATAATCACCAGCAACATGTTTGGTGATATCATCTCCGATGAAGCCTCCGTTATCCCAGGTTCCTTGGGTTTGTTGCCATCTGCGTCCTTGGCCTCTTTGCCAGACAAGAACACCGCATTTGGTTTGTACGAACCATGCCACGGTTCTGCTCCAGATTTGCCAAAGAATAAGGTTGACCCTATCGCCACTATCTTGTCTGCTGCAATGATGTTGAAATTGTCATTGAACTTGCCTGAAGAAGGTAAGGCCATTGAAGATGCAGTTAAAAAGGTTTTGGATGCAGGTATCAGAACTGGTGATTTAGGTGGTTCCAACAGTACCACCGAAGTCGGTGATGCTGTCGCCGAAGAAGTTAAGAAAATCCTTGCTTAAAAAGATTCTCTTTTTTTATGATATTTGTACATAAACTTTATAAATGAAATTCATAATAGAAACGACACGAAATTACAAAATGGAATATGTTCATAGGGTAGACGAAACTATATACGCAATCTACATACATTTATCAAGAAGGAGAAAAAGGAGGATAGTAAAGGAATACAGGTAAGCAAATTGATACTAATGGCTCAACGTGATAAGGAAAAAGAATTGCACTTTAACATTAATATTGACAAGGAGGAGGGCACCACACAAAAAGTTAGGTGTAACAGAAAATCATGAAACTACGATTCCTAATTTGATATTGGAGGATTTTCTCTAAAAAAAAAAAAATACAACAAATAAAAAACACTCAATGACCTGACCATTTGATGGAGTTTAAGTCAATACCTTCTTGAAGCATTTCCCATAATGGTGAAAGTTCCCTCAAGAATTTTACTCTGTCAGAAACGGCCTTACGACGTAGTCGA"),
                id = "leu2",
                name = "leu2",
                description = "leucine 2 auxotroph suppression gene",
                features = [SeqFeature(FeatureLocation(0, 2235), type = "LEU_AUX", strand=1)])

met = SeqRecord(Seq("GCCATCCTCATGAAAACTGTGTAACATAATAACCGAAGTGTCGAAAAGGTGGCACCTTGTCCAATTGAACACGCTCGATGAAAAAAATAAGATATATATAAGGTTAAGTAAAGCGTCTGTTAGAAAGGAAGTTTTTCCTTTTTCTTGCTCTCTTGTCTTTTCATCTACTATTTCCTTCGTGTAATACAGGGTCGTCAGATACATAGATACAATTCTATTACCCCCATCCATACAATGCCATCTCATTTCGATACTGTTCAACTACACGCCGGCCAAGAGAACCCTGGTGACAATGCTCACAGATCCAGAGCTGTACCAATTTACGCCACCACTTCTTATGTTTTCGAAAACTCTAAGCATGGTTCGCAATTGTTTGGTCTAGAAGTTCCAGGTTACGTCTATTCCCGTTTCCAAAACCCAACCAGTAATGTTTTGGAAGAAAGAATTGCTGCTTTAGAAGGTGGTGCTGCTGCTTTGGCTGTTTCCTCCGGTCAAGCCGCTCAAACCCTTGCCATCCAAGGTTTGGCACACACTGGTGACAACATCGTTTCCACTTCTTACTTATACGGTGGTACTTATAACCAGTTCAAAATCTCGTTCAAAAGATTTGGTATCGAGGCTAGATTTGTTGAAGGTGACAATCCAGAAGAATTCGAAAAGGTCTTTGATGAAAGAACCAAGGCTGTTTATTTGGAAACCATTGGTAATCCAAAGTACAATGTTCCGGATTTTGAAAAAATTGTTGCAATTGCTCACAAACACGGTATTCCAGTTGTCGTTGACAACACATTTGGTGCCGGTGGTTACTTCTGTCAGCCAATTAAATACGGTGCTGATATTGTAACACATTCTGCTACCAAATGGATTGGTGGTCATGGTACTACTATCGGTGGTATTATTGTTGACTCTGGTAAGTTCCCATGGAAGGACTACCCAGAAAAGTTCCCTCAATTCTCTCAACCTGCCGAAGGATATCACGGTACTATCTACAATGAAGCCTACGGTAACTTGGCATACATCGTTCATGTTAGAACTGAACTATTAAGAGATTTGGGTCCATTGATGAACCCATTTGCCTCTTTCTTGCTACTACAAGGTGTTGAAACATTATCTTTGAGAGCTGAAAGACACGGTGAAAATGCATTGAAGTTAGCCAAATGGTTAGAACAATCCCCATACGTATCTTGGGTTTCATACCCTGGTTTAGCATCTCATTCTCATCATGAAAATGCTAAGAAGTATCTATCTAACGGTTTCGGTGGTGTCTTATCTTTCGGTGTAAAAGACTTACCAAATGCCGACAAGGAAACTGACCCATTCAAACTTTCTGGTGCTCAAGTTGTTGACAATTTAAAGCTTGCCTCTAACTTGGCCAATGTTGGTGATGCCAAGACCTTAGTCATTGCTCCATACTTCACTACCCACAAACAATTAAATGACAAAGAAAAGTTGGCATCTGGTGTTACCAAGGACTTAATTCGTGTCTCTGTTGGTATCGAATTTATTGATGACATTATTGCAGACTTCCAGCAATCTTTTGAAACTGTTTTCGCTGGCCAAAAACCATGAGTGTGCGTAATGAGTTGTAAAATTATGTATAAACCTACTTTCTCTCACAAG"),
                id = "met15",
                name = "met15",
                description = "methionine auxotroph suppression gene",
                features = [SeqFeature(FeatureLocation(0, 1620), type = "MET_AUX", strand=1)])


################################################################################################
################################################################################################
#
# promoter Sequences
#
#
SP04 = SeqRecord(Seq("TTATTTAGAGTTGACTGTGTTCTCGTCTTCCCATAGCGTGGAAACTGGGGAGTAGCGGA"),
                 name="SP04",
                 description="StrongPromoter04",
                 features = [SeqFeature(FeatureLocation(0, 59), type="promoter", strand=1)])

SP05 = SeqRecord(Seq("CTCCGGAAAGTTGACTTCGAGGGGCCTCATTTATAGAGTGATACCGCCAGAATGCTGGT"),
                 name="SP05",
                 description="StrongPromoter05",
                 features = [SeqFeature(FeatureLocation(0, 59), type="promoter", strand=1)])

SP06 = SeqRecord(Seq("ACTGTTGTAATTGACGAAAGGGGACTTCAAGCGTAGAGTCAACGCTTAAGGCAGGCCAA"),
                name="SP06",
                description="StrongPromoter06",
                features = [SeqFeature(FeatureLocation(0, 59), type="promoter", strand=1)])

SP10 = SeqRecord(Seq("CCCTCGCTGGTTGACACAGTTAGTCAGATTGCCTACGATTTCGTTATCACGCGATTAGC"),
                name="SP10",
                description="StrongPromoter13",
                features = [SeqFeature(FeatureLocation(0, 59), type="promoter", strand=1)])

SP12 = SeqRecord(Seq("TATGCGTTGCTTGACCAAACCTATGTATAGGGATAGGGTTGGTCGGCCCAGTAGTTATC"),
                name="SP12",
                description="StrongPromoter12",
                features = [SeqFeature(FeatureLocation(0, 59), type="promoter", strand=1)])

SP13 = SeqRecord(Seq("CGATGTCGCGTTGACAACGCTGAGGTGAGGGGTTAGGGTAGGTGATCACGGAAGTAGCC"),
                name="SP13",
                description="StrongPromoter13",
                features = [SeqFeature(FeatureLocation(0, 59), type="promoter", strand=1)])

SP14 = SeqRecord(Seq("GTGGGTCGTATTGACGATGTATTTGGGGAGTGTTAGGATATGCGCAATAAAGATGTCGT"),
                name="SP14",
                description="StrongPromoter14",
                features = [SeqFeature(FeatureLocation(0, 59), type="promoter", strand=1)])

SP15 = SeqRecord(Seq("AATTGCCCACTTGACGTTGAGAGTGAAGCAATATAGGTTAACCTGCCACGTGTTCGTTA"),
                name="SP15",
                description="StrongPromoter15",
                features = [SeqFeature(FeatureLocation(0, 59), type="promoter", strand=1)])

SP16 = SeqRecord(Seq("GGTGATCAGTTTGACGCCCCTGACATCGGGTGGTAGCGTCGCTGGTGGTCCGTCAAGCC"),
                name="SP16",
                description="StrongPromoter16",
                features = [SeqFeature(FeatureLocation(0, 59), type="promoter", strand=1)])

SP17 = SeqRecord(Seq("TCGGTGGACTTTGACGTACGGGGTCCCGGTCTGTAGCATGAGGACCGAATACGGGTCCT"),
                name="SP17",
                description="StrongPromoter17",
                features = [SeqFeature(FeatureLocation(0, 59), type="promoter", strand=1)])

SP18 = SeqRecord(Seq("GGCGTAGGGTTTGACTAAAATAATATTCCATTATAGAGTCGTCAGGTCTGGTCTACACA"),
                name="SP18",
                description="StrongPromoter18",
                features = [SeqFeature(FeatureLocation(0, 59), type="promoter", strand=1)])

SP19 = SeqRecord(Seq("TACATCGCAGTTGACAAGCCGGCTTTCTGATCATACGGTGTGACGAAGTTGCTGGGCGT"),
                name="SP19",
                description="StrongPromoter19",
                features = [SeqFeature(FeatureLocation(0, 59), type="promoter", strand=1)])

SP21 = SeqRecord(Seq("GTCTGTGTCCTTGACGCCGCGGGGGTGTAGAGATAGAGTCTGGCCTGCTCTGGCCAACC"),
                name="SP21",
                description="StrongPromoter21",
                features = [SeqFeature(FeatureLocation(0, 59), type="promoter", strand=1)])

SP23 = SeqRecord(Seq("TCACATGTGTTTGACCCCGTGGGTCTGAGTTGCTACCATCAAAGCTGCTCTGGCCAGCC"),
                name="SP23",
                description="StrongPromoter23",
                features = [SeqFeature(FeatureLocation(0, 59), type="promoter", strand=1)])

################################################################################################
################################################################################################
#
# Ribosomal Binding Sites
#
#
RBS04 = SeqRecord(Seq("CTAGGAGGGTTGTGA"),
                name="RBS04",
                description="Ribosomal Binding Site 04",
                features = [SeqFeature(FeatureLocation(0, 15), type="RBS", strand=1)])

RBS05 = SeqRecord(Seq("TCAGGAGGTGAATTA"),
                name="RBS05",
                description="Ribosomal Binding Site 05",
                features = [SeqFeature(FeatureLocation(0, 15), type="RBS", strand=1)])

RBS06 = SeqRecord(Seq("TCAGGAGGTAACACA"),
                name="RBS06",
                description="Ribosomal Binding Site 06",
                features = [SeqFeature(FeatureLocation(0, 15), type="RBS", strand=1)])

RBS07 = SeqRecord(Seq("TGAGGAGGACCCTTT"),
                name="RBS07",
                description="Ribosomal Binding Site 07",
                features = [SeqFeature(FeatureLocation(0, 15), type="RBS", strand=1)])

RBS08 = SeqRecord(Seq("ACAGGAGGGGCGATC"),
                name="RBS08",
                description="Ribosomal Binding Site 08",
                features = [SeqFeature(FeatureLocation(0, 15), type="RBS", strand=1)])

RBS09 = SeqRecord(Seq("TGAGGAGGAGCGATG"),
                name="RBS09",
                description="Ribosomal Binding Site 09",
                features = [SeqFeature(FeatureLocation(0, 15), type="RBS", strand=1)])

RBS10 = SeqRecord(Seq("TAAGGAGGTAAAGGC"),
                name="RBS10",
                description="Ribosomal Binding Site 10",
                features = [SeqFeature(FeatureLocation(0, 15), type="RBS", strand=1)])

RBS11 = SeqRecord(Seq("GCAGGAGGACTGTAT"),
                name="RBS11",
                description="Ribosomal Binding Site 11",
                features = [SeqFeature(FeatureLocation(0, 15), type="RBS", strand=1)])

RBS12 = SeqRecord(Seq("GCAGGAGGACTGTAT"),
                name="RBS12",
                description="Ribosomal Binding Site 12",
                features = [SeqFeature(FeatureLocation(0, 15), type="RBS", strand=1)])

RBS13 = SeqRecord(Seq("TGAGGAGGGTTGAAT"),
                name="RBS13",
                description="Ribosomal Binding Site 13",
                features = [SeqFeature(FeatureLocation(0, 15), type="RBS", strand=1)])

RBS14 = SeqRecord(Seq("AAAGGAGGGATGGGT"),
                name="RBS14",
                description="Ribosomal Binding Site 14",
                features = [SeqFeature(FeatureLocation(0, 15), type="RBS", strand=1)])

RBS15 = SeqRecord(Seq("CAAGGAGGTTCCAAT"),
                name="RBS15",
                description="Ribosomal Binding Site 15",
                features = [SeqFeature(FeatureLocation(0, 15), type="RBS", strand=1)])

RBS16 = SeqRecord(Seq("GTAGGAGGTTTGGTA"),
                name="RBS16",
                description="Ribosomal Binding Site 16",
                features = [SeqFeature(FeatureLocation(0, 15), type="RBS", strand=1)])

RBS17 = SeqRecord(Seq("TGAGGAGGCGATAGA"),
                name="RBS17",
                description="Ribosomal Binding Site 17",
                features = [SeqFeature(FeatureLocation(0, 15), type="RBS", strand=1)])

RBS18 = SeqRecord(Seq("CAAGGAGGACGTCTT"),
                name="RBS18",
                description="Ribosomal Binding Site 18",
                features = [SeqFeature(FeatureLocation(0, 15), type="RBS", strand=1)])


################################################################################################
################################################################################################
#
# Insulator Sequences
# names after the SP sequence they are paired with
#
insSP04 = SeqRecord(Seq("GGGACAACATTCTAACACCCAAGAGTTGCTTCGCGCACAGCCAGCCTTGG"),
                name="insSP04",
                description="Insulator Sequence SP04",
                features = [SeqFeature(FeatureLocation(0, 50), type="stuffer", strand=1)])

insSP05 = SeqRecord(Seq("CATCTCCCCGACGGTCACAATGGGGCGCTCCTAGGCTCTTCCGCTTGGCC"),
                name="insSP05",
                description="Insulator Sequence SP05",
                features = [SeqFeature(FeatureLocation(0, 50), type="stuffer", strand=1)])

insSP06 = SeqRecord(Seq("GTGACATCAATAGAGGGGCAGGCTCTCGCTATGTCTCTCAAGGGCATGCC"),
                name="insSP06",
                description="Insulator Sequence SP06",
                features = [SeqFeature(FeatureLocation(0, 50), type="stuffer", strand=1)])

insSP10 = SeqRecord(Seq("ACCAAGGTGGGATAGCCCCGTAGCCGACTACGCCTCCCTCTGGGCCTTAC"),
                name="insSP10",
                description="Insulator Sequence SP10",
                features = [SeqFeature(FeatureLocation(0, 50), type="stuffer", strand=1)])

insSP12 = SeqRecord(Seq("GTCTCTTTGCCGACTAATGCGAACAACCACACCATAGCGATTCGTCGGGG"),
                name="insSP12",
                description="Insulator Sequence SP12",
                features = [SeqFeature(FeatureLocation(0, 50), type="stuffer", strand=1)])

insSP13 = SeqRecord(Seq("ACCCGGACGCGTGGCACCTCAGGAGGCGGCCGGAGGGGGGATGTCCTCTG"),
                name="insSP13",
                description="Insulator Sequence SP13",
                features = [SeqFeature(FeatureLocation(0, 50), type="stuffer", strand=1)])

insSP14 = SeqRecord(Seq("CGCTGATGGCGCTCACCAGGGCAGTTTTATCCAGCCCCATAGGTGTTCAC"),
                name="insSP14",
                description="Insulator Sequence SP14",
                features = [SeqFeature(FeatureLocation(0, 50), type="stuffer", strand=1)])

insSP15 = SeqRecord(Seq("GGAGGGTGCCCTCCCAGCACATGAACTGGGGCCACAGAGACGCGTGTAGA"),
                name="insSP15",
                description="Insulator Sequence SP15",
                features = [SeqFeature(FeatureLocation(0, 50), type="stuffer", strand=1)])

insSP16 = SeqRecord(Seq("GGCAGGGGTAGGACCATCGGTAGTAGGGATAGTGCGGAAGCTCGCTGACC"),
                name="insSP16",
                description="Insulator Sequence SP16",
                features = [SeqFeature(FeatureLocation(0, 50), type="stuffer", strand=1)])

insSP17 = SeqRecord(Seq("CCCCAGCAGTTCGTCTCGCGTGCGGGCCTTCGCTACCTGCACGATACCTG"),
                name="insSP17",
                description="Insulator Sequence SP17",
                features = [SeqFeature(FeatureLocation(0, 50), type="stuffer", strand=1)])

insSP18 = SeqRecord(Seq("AACCGCGTCTCCACGACCGGCGCTCGATTCAACTTCGCCGACGTGACGAC"),
                name="insSP18",
                description="Insulator Sequence SP18",
                features = [SeqFeature(FeatureLocation(0, 50), type="stuffer", strand=1)])

insSP19 = SeqRecord(Seq("AGATGTAGCGGCTAGGAACGCAACAAGCATCGACGAACGGCCTCGAATAG"),
                name="insSP19",
                description="Insulator Sequence SP19",
                features = [SeqFeature(FeatureLocation(0, 50), type="stuffer", strand=1)])

insSP21 = SeqRecord(Seq("GGCTTTCACTGCGAGGTGGCTAAACGAGTACACGCTCGTTACTTGTCATA"),
                name="insSP21",
                description="Insulator Sequence SP21",
                features = [SeqFeature(FeatureLocation(0, 50), type="stuffer", strand=1)])

insSP23 = SeqRecord(Seq("GATACTGCCATGGACGACTACCCATCCCTCTGGGCCTCAGACAGCCGGAT"),
                name="insSP23",
                description="Insulator Sequence SP23",
                features = [SeqFeature(FeatureLocation(0, 50), type="stuffer", strand=1)])


################################################################################################
################################################################################################
#
# COmbined Promoter Sequences
#
#
#
# Dual Promoter Sequences:
#SNP11 = rbs104 (r) + SP15 (r) + LEU + Sp16 + rbs105
#SNP12 = rbs107 (r) + SP14 (r) + insSP14 (r) + insSP13 + Sp13 + rbs108
#SNP13 = rbs109 (r) + SP19 (r) + MET + Sp12 + rbs110
#SNP14 = rbs106 (r) + SP21 (r) + insSP21 (r) + insSP23 + Sp23 + rbs112
#SNP15 = rbs113 (r) + SP17 (r) + insSP17 (r) + insSP18 + Sp18 + rbs116
#SNP16 = rbs117 (r) + SP06 (r) + insSP06 (r) + insSP10 + Sp10 + rbs118
#SNP17 = rbs114 (r) + SP04 (r) + insSP04 (r) + insSP05 + Sp05 + rbs115
SNP11 = RBS04.reverse_complement() + SP15.reverse_complement() + leu + SP16 + RBS05
SNP12 = RBS07.reverse_complement() + SP14.reverse_complement() + insSP14.reverse_complement() + insSP13 + SP13 + RBS08
SNP13 = RBS09.reverse_complement() + SP19.reverse_complement() + met + SP12 + RBS10
SNP14 = RBS06.reverse_complement() + SP21.reverse_complement() + insSP21.reverse_complement() + insSP23 + SP23 + RBS12
SNP15 = RBS13.reverse_complement() + SP17.reverse_complement() + insSP17.reverse_complement() + insSP18 + SP18 + RBS16
SNP16 = RBS17.reverse_complement() + SP06.reverse_complement() + insSP06.reverse_complement() + insSP10 + SP10 + RBS18
SNP17 = RBS14.reverse_complement() + SP04.reverse_complement() + insSP04.reverse_complement() + insSP05 + SP05 + RBS15

#set auxotrophic information:
SNP11.id = "SNP11"
SNP12.id = "SNP12"
SNP13.id = "SNP13"
SNP14.id = "SNP14"
SNP15.id = "SNP15"
SNP16.id = "SNP16"
SNP17.id = "SNP17"


# in order of effectiveness as auxotrophic markers
selective_promoters = [SNP11, SNP13]
nonselective_promoters = [SNP12, SNP14, SNP15, SNP16, SNP17]


################################################################################################
################################################################################################
#
# CRISPR-Related Sequences
#
#
# CRISPR Direct Repeat
crispr_DR = SeqRecord(Seq("GTTTTAGAGCTATGCTGTTTTGAATGGTCCCAAAAC"),
                name="CRISPR_Direct_Repeat",
                description="CRISPR_Direct_Repeat",
                features = [SeqFeature(FeatureLocation(0, 36), type="CRISPR_Direct_Repeat", strand=1)])

gateway_5prime = SeqRecord(Seq("CTTTGGTCTCACCAAAAC"),
                name="Gateway 5prime",
                description="Gateway 5prime",
                # features = [SeqFeature(FeatureLocation(5, 10), type="BSAI_Recognition", strand=1),
                #             SeqFeature(FeatureLocation(12, 18), type="BSAI_Cut_Site", strand=1)])
                features = [SeqFeature(FeatureLocation(4, 10), type="BSAI_Recognition", strand=1),
                            SeqFeature(FeatureLocation(11, 18), type="BSAI_Cut_Site", strand=1)])

gateway_3prime = SeqRecord(Seq("GTTTTAGAGAGAGACCTTTC"),
                name="Gateway 3prime",
                description="Gateway 3prime",
                # features = [SeqFeature(FeatureLocation(1, 9), type="BSAI_Cut_Site", strand=1),
                #             SeqFeature(FeatureLocation(11, 17), type="BSAI_Recognition", strand=1)])
                features = [SeqFeature(FeatureLocation(0, 9), type="BSAI_Cut_Site", strand=1),
                            SeqFeature(FeatureLocation(10, 16), type="BSAI_Recognition", strand=1)])
