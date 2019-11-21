from collections import defaultdict
import collections
from scipy import stats
from scipy.stats import mannwhitneyu
import numpy
import matplotlib
import matplotlib.pyplot as plt
import pylab
from itertools import islice
from operator import itemgetter



def average_standard_dev(lengths):
    """function to return the avaerage and stadard deviation
    for a list of number.
    Uses Numpy to do the calculations.
    Take in a list of lens
    returns the mean and standard deviation of the list"""
    the_mean = sum(lengths) / float(len(lengths))
    standard_dev = numpy.std(lengths)
    return the_mean, standard_dev


def stats_on_list_of_sizes(db_lens, assemb_lens):
    """function to perform stats on two lists of seq lens.
    Returns as a tab separeated string:
    as_skew,
    db_skew,
    ttest,
    Man_u_value,
    Man_p_value"""
    as_skew = ('normal skewtest assemb_lens = %6.3f pvalue = %6.4f' %
               stats.skewtest(assemb_lens))
    db_skew = ('normal skewtest db_lens = %6.3f pvalue = %6.4f' %
               stats.skewtest(db_lens))
    ttest = ('t-statistic = %6.3f pvalue = %6.4f' %
             stats.ttest_ind(db_lens, assemb_lens))
    Man_u_value, Man_p_value = mannwhitneyu(db_lens, assemb_lens,
                                            alternative="two-sided")
    outdata = "\t".join([as_skew,
                         db_skew,
                         ttest,
                         str(Man_u_value),
                         str(Man_p_value)])
    return outdata


def plot_seq_len_histograms(db, assembled, outname):
    """takes in a 2 lists of number and plots a histogram.
    Plots two histograms on a graph """
    # the histogram of the data
    fig = plt.figure(figsize=(10, 8), dpi=1200)
    ax1 = fig.add_subplot(1, 2, 1)  # 1x2 grid, position 1
    ax2 = fig.add_subplot(1, 2, 2)  # 1x2 grid, position 2
    # graph1
    rects1 = ax1.hist(db, facecolor='green', alpha=0.6)
    ax1.set_xlabel('Database sequence lengths')
    ax1.set_ylabel('Count')
    ax1.grid(True)
    ax1.set_title("Histogram of database sequence lengths")
    # graph 2
    rects2 = ax2.hist(assembled, facecolor='green', alpha=0.6)
    ax2.set_xlabel('Assembled sequence lengths')
    ax2.set_ylabel('Count')
    ax2.grid(True)
    ax2.set_title("Histogram of assembled sequence lengths")
    fig.tight_layout()
    fig
    outpng = os.path.join(outname)
    pylab.savefig(outpng)
    pylab.close()


def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))

def open_file_retunrPI(infile):
    """def to open the pi ouput from vcf tools
    return a list of the vals and other dictionaries. """
    data_list = []
    scaff_bin_to_PI = defaultdict(int)
    scaff_bin_to_varients = defaultdict(int)
    with open(infile) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if line.startswith("CHROM"):
                continue
            if not line.strip():
                    continue  #  if the last line is blank
            CHROM, BIN_START, BIN_END, N_VARIANTS, PI = line.rstrip().split()
            PI = float(PI)
            scaff_bin = CHROM + "_" + BIN_START
            data_list.append(PI)
            scaff_bin_to_PI[scaff_bin] = PI
    return data_list, scaff_bin_to_PI, scaff_bin_to_varients
      

    
# lets parse some data
cress_pi_data, cress_scaff_bin_to_PI, \
               cress_scaff_bin_to_varients = open_file_retunrPI("M.cerasi_cress.10000.sitepi.windowed.pi")
cherry_pi_data, cherry_scaff_bin_to_PI,\
               cherry_scaff_bin_to_varients = open_file_retunrPI("M.cerasi_cherry.10000.sitepi.windowed.pi")
cleavers_pi_data, cleavers_scaff_bin_to_PI, \
                  cleavers_scaff_bin_to_varients = open_file_retunrPI("M.cerasi_cleavers.10000.sitepi.windowed.pi")

# avs and sd of the data:
av, sd = average_standard_dev(cress_pi_data)

print("cress av and sd = ", av, sd)

av, sd = average_standard_dev(cherry_pi_data)

print("cherry av and sd = ", av, sd)

av, sd = average_standard_dev(cleavers_pi_data)

print("cleavers av and sd = ", av, sd)
# get some stats on the data
cherry_vs_creess_stats = stats_on_list_of_sizes(cherry_pi_data, cress_pi_data)
print("cherry_vs_creess_stats: ", cherry_vs_creess_stats)

# nope not enough memory
#plot_seq_len_histograms(cherry_pi_data, cress_pi_data, "cress_cherry.png")

cherry_vs_cleav_stats = stats_on_list_of_sizes(cherry_pi_data, cleavers_pi_data)
print("cherry_vs_cleav_stats:", cherry_vs_cleav_stats)

cress_vs_cleav_stats = stats_on_list_of_sizes(cress_pi_data, cleavers_pi_data)
print("cress_vs_cleav_stats:", cress_vs_cleav_stats)

data = [cherry_pi_data, cleavers_pi_data, cress_pi_data]
fig7, ax7 = plt.subplots()
matplotlib.pyplot.yscale("log")
ax7.set_title('Boxplot of PI ')
ax7.boxplot(data,  notch=True)

#plt.show()

# compare SNP hetero / homo  rations
cress_hetero_ratio = [0.45,0.62,0.45]
cleaver_hetero_ratio = [0.64,0.48,0.45]
cherry_hetero_ratio = [0.23,0.26,0.26]

ttest = ('t-statistic = %6.3f pvalue = %6.4f' %
             stats.ttest_ind(cress_hetero_ratio, cleaver_hetero_ratio))
print("gal_vs_creess_Hetero_ratio: ", ttest)


ttest = ('t-statistic = %6.3f pvalue = %6.4f' %
             stats.ttest_ind(cress_hetero_ratio, cherry_hetero_ratio))
print("cherry_vs_creess_Hetero_ratio: ", ttest)

ttest = ('t-statistic = %6.3f pvalue = %6.4f' %
             stats.ttest_ind(cleaver_hetero_ratio, cherry_hetero_ratio))
print("cherry_vs_cleavers_Hetero_ratio: ", ttest)

# Cress: find the most SNP heavy region:
cress_scaff_bin_to_PI2 = collections.OrderedDict(sorted(cress_scaff_bin_to_PI.items(),
                                                             key=itemgetter(1), reverse=False))
cress_scaff_bin_to_varients2 = collections.OrderedDict(sorted(cress_scaff_bin_to_varients.items(),
                                                             key=itemgetter(1), reverse=False))

n_items = take(3, cress_scaff_bin_to_PI2.items())

print("cress, hiest PI:", n_items)
n_items = take(10, cress_scaff_bin_to_varients2.items())
print("cress, most varients:", n_items)

# cherry: find the most SNP heavy region:
cherry_scaff_bin_to_PI2 = collections.OrderedDict(sorted(cherry_scaff_bin_to_PI.items(),
                                                             key=itemgetter(1), reverse=False))
cherry_scaff_bin_to_varients2 = collections.OrderedDict(sorted(cherry_scaff_bin_to_varients.items(),
                                                             key=itemgetter(1), reverse=False))

n_items = take(10, cherry_scaff_bin_to_PI2.items())

print("cherry, hiest PI:", n_items)
n_items = take(10, cherry_scaff_bin_to_varients2.items())
print("cherry, most varients:", n_items)

# cleavers: find the most SNP heavy region:
cleavers = collections.OrderedDict(sorted(cress_scaff_bin_to_PI.items(),
                                                             key=itemgetter(1), reverse=False))
cress_scaff_bin_to_varients2 = collections.OrderedDict(sorted(cress_scaff_bin_to_varients.items(),
                                                             key=itemgetter(1), reverse=False))

n_items = take(3, cleavers.items())

print("cress, hiest PI:", n_items)
n_items = take(3, cress_scaff_bin_to_varients2.items())
print("cress, most varients:", n_items)
