import os
import re

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib



plt.rcParams["figure.facecolor"] = "white"
plt.rcParams["axes.facecolor"] = "white"
plt.rcParams["savefig.facecolor"] = "white"

plt.rcParams.update({
	"text.usetex" : True,
    "font.family" : "serif",
    "font.serif" : ["Computer Modern Serif"],
})

def combine_down_histograms_snr(background_dir):
    """
    This will reload and combine down the histograms
    """

    # get a list of all the files in tbe directory containing
    # the histograms
    list_all_files = os.listdir(background_dir)
    # one of these is the information on the bin edges
    bin_edges_fname = [i for i in list_all_files if 'bin_edges.txt' in i]


    # the rest are actual histograms
    histogram_fnames = [i for i in list_all_files if 'zerolags' in i]

    bin_edges = np.loadtxt(background_dir+'/'+bin_edges_fname[0])

    total_histogram = [0]*(len(bin_edges)-1)

    for file in histogram_fnames:
        this_histogram = np.loadtxt(background_dir+'/'+file)
        for this_bin in range(len(bin_edges)-1):
            total_histogram[this_bin] += this_histogram[this_bin]

    return total_histogram, list(bin_edges)

def combine_down_histograms_chisq(background_dir):
    """
    This will reload and combine down the histograms

    """

    # get a list of all the files in tbe directory containing
    # the histograms
    list_all_files = os.listdir(background_dir)
    # one of these is the information on the bin edges
    bin_edges_fname = [i for i in list_all_files if 'bin_edges_chisq.txt' in i]


    # the rest are actual histograms
    histogram_fnames = [i for i in list_all_files if 'chisq_' in i]
    #print(histogram_fnames)

    bin_edges = np.loadtxt(background_dir+'/'+bin_edges_fname[0])

    total_histogram = [0]*(len(bin_edges)-1)

    for file in histogram_fnames:
        this_histogram = np.loadtxt(background_dir+'/'+file)
        for this_bin in range(len(bin_edges)-1):
            total_histogram[this_bin] += this_histogram[this_bin]
    #print(len(total_histogram))

    return total_histogram, list(bin_edges)


#################
# Example usage #
#################

# # for the SNR histogram of all the background
# # suppose your background files all contain the string 'new_background':
# BACKGROUND_DIRS = [f for f in os.listdir(CWD) if 'new_background' in f]

# hist_snr = [0]*100 #100 is the number of bins, can change
# for this_file in BACKGROUND_DIRS:
#     this_hist, bin_edges = combine_down_histograms_snr(CWD+'/'+this_file)
#     hist_snr = [sum(x) for x in zip(hist_snr_wf1, this_hist)]

# hist_chisq_wf1 = [0]*200 #200 is the number of bins, can change
# for this_file in BACKGROUND_DIRS:
#     this_hist, bin_edges_chisq = combine_down_histograms_chisq(CWD+'/'+this_file)
#     hist_chisq = [sum(x) for x in zip(hist_chisq, this_hist)]

# # then to plot, this just does the chisq, but you can also change to SNR appropriately:
# chisq_bins_unlogged = [10**ii for ii in bin_edges_chisq]
# snr_bins_unlogged = [10**ii for ii in bin_edges]

# plt.figure(figsize = (9, 4))
# plt.hist(chisq_bins_unlogged[:-1], chisq_bins_unlogged, 
#         weights = hist_chisq, color = 'DeepSkyBlue',
#         histtype = 'step', lw = 2, 
#         label = 'Chisq distribution of background triggers')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlim([1E-1, 1E2])
# plt.show()
