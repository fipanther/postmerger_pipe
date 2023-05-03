import data_preprocess as pre
import trigger_finder as pks
import pycbc.psd
import pycbc.filter
import pycbc.waveform
from pycbc.vetoes import power_chisq
from pycbc.types import TimeSeries, FrequencySeries

from scipy.signal.windows import tukey

import numpy as np
import matplotlib.pyplot as plt

import os


#############
# WAVEFORMS #
#############

# note that from here to line 133 should really 
# go into a separate module/YAML file combo

####
# Template 1 parameters
####
# # here we will generate the template - this has an overlap of ~65%
wf_order = 3
# the template is going to be renormalized regardless
log_amplitude = 1 
t_0 = 0
tukey_rolloff = 0#0.01
w_0 = 0.35
w_1 = 0.56
f = [3305, 2540, 1870]
T = [3.0, 1.25, 1.5] # note that T should be specified in ms
phi = [1.74, 0.68, -0.2]


####
# Template 2
####
# here we will generate the template - this has an overlap of ~65% with Template 1
# wf_order = 3
# # the template is going to be renormalized regardless
# log_amplitude = 1 
# t_0 = 0
# tukey_rolloff = 0
# w_0 = 0.4
# w_1 = 0.54
# f = [3250, 2500, 1890]
# T = [3.5, 1.9, 0.1] # note that T should be specified in ms
# phi = [1.57, 0.7, 1.57]

###########################################################
#  Mostly PyCBC boiler-plate, with some heavy inspiration #
#  from Paul Easter's post-merger PE code                 #
###########################################################

def freq_wf(frequency, damping, central_freq, phase, amplitude, weight):
    """
    returns the h+ component of the fd waveform of a single damped sinusoid
    """
    
    # plus component
    omega_0 = 2 * np.pi * central_freq
    damp = 1 / (damping/1000.0)
    A = np.exp(1.0j * phase)/(damp - (1.0j * (omega_0 - (2*np.pi*frequency))))
    B = np.exp(-1.0j * phase)/(damp + (1.0j * (omega_0 +  (2*np.pi*frequency))))
    
    hplus = amplitude * weight * 0.5 * (A + B)
    hcross = amplitude * weight * (0.5/1.0j) * (A - B)
    
    return hplus, hcross

def f_domain_n_component_damped_sinusoid(frequency, **args):

    amplitude = args['log_amplitude']
    hplus = np.zeros(len(frequency), dtype=np.complex_)
    hcross = np.zeros(len(frequency), dtype=np.complex_)
    
    for waveform_number in range(3):
        if waveform_number == 2:
            weight = (1 - args['w_0'] - args['w_1'])
        else:
            weight = args[f'w_{waveform_number}']
        A = 10 ** amplitude
        f_0 = args[f'f_{waveform_number}']
        T = args[f'T_{waveform_number}']
        phi = args[f'phi_{waveform_number}']
        hplus += freq_wf(frequency, T, f_0, phi, A, weight)[0]
        hcross += freq_wf(frequency, T, f_0, phi, A, weight)[1]
    
    return hplus, hcross

def waveform_3_comp_fd(frequency_in = None, **args):
    # waveform model from Panther+2023, in the frequency domain
    delta_f = args['delta_f']
    frequency = np.linspace(min(frequency_in), max(frequency_in), len(frequency_in))

    A = args['log_amplitude']

    w_0 = args['w_0']
    w_1 = args['w_1']

    f_0 = args['f_0']
    f_1 = args['f_1']
    f_2 = args['f_2']

    T_0 = args['T_0']
    T_1 = args['T_1']
    T_2 = args['T_2']

    phi_0 = args['phi_0']
    phi_1 = args['phi_1']
    phi_2 = args['phi_2']

    hplus, hcross = f_domain_n_component_damped_sinusoid(frequency, log_amplitude = A, 
                                                        w_0 = w_0, w_1 = w_1,
                                                        f_0 = f_0, f_1 = f_1, f_2 = f_2, 
                                                        T_0 = T_0, T_1 = T_1, T_2 = T_2, 
                                                        phi_0 = phi_0, phi_1 = phi_1, phi_2 = phi_2)                               
    
    #waveform = hplus + hcross

    hplus_out = FrequencySeries(hplus, delta_f= delta_f)
    hcross_out = FrequencySeries(hcross, delta_f= delta_f)
    
    return hplus_out, hcross_out

pycbc.waveform.add_custom_waveform('waveform_3_comp_fd', waveform_3_comp_fd, 'frequency', force=True)

###################
#   USER INPUTS   #
###################

#------- where do you want to store your outputs? -------#
# where do you want to put the results
CWD = os.getcwd()
background_dirname = CWD+'/new_background_fd_template_wf1_chunk1_chisq/'

#------- where did you store the strain data? ---------#
DATA_DIR = '/Users/fipanther/Documents/Work/postmerger/ligo_frames/'
filename_L1 ='L1/L-L1_GWOSC_O3a_16KHZ_R1-1238175744-4096.hdf5'
filename_H1 ='H1/H-H1_GWOSC_O3a_16KHZ_R1-1238175744-4096.hdf5'

# how many 128s chunks do you wait to drop a file?
nchunks = 10

# what do you consider to be the threshold for the 'loudest' triggers you might need to investigate?
loudest_thresh = 6.
loudest_thresh_fname = 'significant_events.txt'


#################
#  START HERE   #
#################

# histogramming setup
# bins are defined in log10(rho_c) and log10(chisq)
# you can change these if you want as you might lose
# some of the loudest triggers
hist_bins = np.linspace(0.25, 2.0, 100)
hist_bins_chisq = np.linspace(-2, 2, 200)



# import the strain from the file and make the 128s chunks
# NOTE: I have not considered that there may be data gaps
# to account for data gaps, we need to use the PyCBC 'gate' function
# and apply DQ vetos
strain_full_L1, t_start_L1, t_end_L1, duration_L1, ts_L1, sample_rate_L1 = pre.import_data_from_file(DATA_DIR, filename_L1, type = 'hdf5')
strain_as_timeseries_L1 = strain_full_L1
chunks_L1 = pre.batch_strain(strain_as_timeseries_L1, t_start_L1, ts_L1)

strain_full_H1, t_start_H1, t_end_H1, duration_H1, ts_H1, sample_rate_H1 = pre.import_data_from_file(DATA_DIR, filename_H1, type = 'hdf5')
strain_as_timeseries_H1 = strain_full_H1


chunks_H1 = pre.batch_strain(strain_as_timeseries_H1, t_start_H1, ts_H1)


# create some empty lists to store your output
triggers_saved = []
idx_saved = []
indiv_saved = []
chisq_saved = []
cmb_chisq_saved = []
idx_over_thresh = []
chunk_over_thresh = []
buffer_no_over_thresh = []
indiv_over_thresh = []
chisq_over_thresh = []

# here we will loop through the chunks

for chunk_i in range(0, min(len(chunks_L1),len(chunks_H1))):
	print('Analysing chunk {}'.format(chunk_i))
	# # precondition each chunk
	# # note that each preconditioned chunk loses 2s of the beginning and end
	# # to remove spectral leakage. 
	# # Thus, we need to keep careful track of the live time
	chunk_1_conditioned_L1 = pre.data_precondition(chunks_L1[chunk_i], ts_L1)
	chunk_1_conditioned_H1 = pre.data_precondition(chunks_H1[chunk_i], ts_H1)

	# We will perform frequency domain filtering here
	stilde_L1 = pycbc.filter.make_frequency_series(chunk_1_conditioned_L1)
	stilde_H1 = pycbc.filter.make_frequency_series(chunk_1_conditioned_H1)

	# #We'll choose 4 seconds PSD samples that are overlapped 50 %
	# # By default pycbc.psd.welch uses the median which is appropriately
	# # biased. This function also windows the data using a Hann window

	# here I assume that ts_L1 == ts_H1 (sample rate) within machine precision
	seg_len = int(4 / ts_L1)
	seg_stride = int(seg_len / 2)

	estimated_psd_L1 = pycbc.psd.welch(chunk_1_conditioned_L1,
	                      seg_len=seg_len,
	                      seg_stride=seg_stride)



	estimated_psd_H1 = pycbc.psd.welch(chunk_1_conditioned_H1,
	                      seg_len=seg_len,
	                      seg_stride=seg_stride)



	# interpolate to the correct delta_f
	psd_interpolated_L1 = pycbc.psd.interpolate(estimated_psd_L1, stilde_L1.delta_f)
	psd_interpolated_H1 = pycbc.psd.interpolate(estimated_psd_H1, chunk_1_conditioned_H1.delta_f)


	# FOR THE FD CODE
	filter_L1, _ = pycbc.waveform.get_fd_waveform(approximant="waveform_3_comp_fd",
                                    frequency_in = stilde_L1.sample_frequencies,
                                    delta_f = stilde_L1.delta_f,
                                    log_amplitude = log_amplitude,
                                    w_0 = w_0 , w_1 = w_1,
                                    f_0 = f[0], f_1 = f[1], f_2 = f[2],
                                    T_0 = T[0], T_1 = T[1], T_2 = T[2],
                                    phi_0 = phi[0], phi_1 = phi[1], phi_2 = phi[2], 
                                    f_lower = 20)

	filter_H1, _ = pycbc.waveform.get_fd_waveform(approximant="waveform_3_comp_fd",
                                    frequency_in = stilde_H1.sample_frequencies,
                                    delta_f = stilde_H1.delta_f,
                                    log_amplitude = log_amplitude,
                                    w_0 = w_0 , w_1 = w_1,
                                    f_0 = f[0], f_1 = f[1], f_2 = f[2],
                                    T_0 = T[0], T_1 = T[1], T_2 = T[2],
                                    phi_0 = phi[0], phi_1 = phi[1], phi_2 = phi[2], 
                                    f_lower = 20)
	
	# then you can actually generate the SNR timeseries
	snr_full_L1 = pycbc.filter.matched_filter(filter_L1, stilde_L1,
	                     psd=psd_interpolated_L1, low_frequency_cutoff=1000, high_frequency_cutoff=4096)
	snr_full_H1 = pycbc.filter.matched_filter(filter_H1, stilde_H1,
	                     psd=psd_interpolated_H1, low_frequency_cutoff=1000, high_frequency_cutoff=4096)

	#  remove the artefacts from the filter with the tukey window 

	snr_full_L1 = snr_full_L1 * tukey(len(snr_full_L1), 0.5 / len(snr_full_L1)/snr_full_L1.sample_rate)	
	snr_full_H1 = snr_full_H1 * tukey(len(snr_full_H1), 0.5 / len(snr_full_L1)/snr_full_L1.sample_rate)

	# once you have it you can start finding triggers

	# EDIT 29/11: add chisq to check vetos
	nbins = 10 # this parameter is set empirically - you will need to tune it 
	chisq_L1 = power_chisq(filter_L1, stilde_L1, nbins, psd_interpolated_L1,
                low_frequency_cutoff=1000,
                high_frequency_cutoff=4096,
                return_bins=False)
	chisq_H1 = power_chisq(filter_H1, stilde_H1, nbins, psd_interpolated_H1,
                low_frequency_cutoff=1000,
                high_frequency_cutoff=4096,
                return_bins=False)

	chisq_L1 /=  (nbins * 2) - 2 # we only care about reduced chisq
	chisq_H1 /=  (nbins * 2) - 2
	# end edit

	search_buffer_list_L1, search_buffer_list_chisq_L1, t_ends_L1 = pks.create_search_buffers(snr_full_L1, chisq_L1, 1, 0)
	search_buffer_list_H1, search_buffer_list_chisq_H1, t_ends_H1 = pks.create_search_buffers(snr_full_H1, chisq_H1, 1, 0)



	for buffer_i in range(0,len(search_buffer_list_L1)):


		# find the best trigger in each second for each IFO
		triggers_L1, trigger_idx_L1, indiv_L1, indiv_chisq_L1 = pks.find_coincident_peaks(search_buffer_list_L1[buffer_i], search_buffer_list_H1[buffer_i], search_buffer_list_chisq_L1[buffer_i], search_buffer_list_chisq_H1[buffer_i], 4, 0, 0.01)
		triggers_H1, trigger_idx_H1, indiv_H1, indiv_chisq_H1 = pks.find_coincident_peaks(search_buffer_list_H1[buffer_i], search_buffer_list_L1[buffer_i], search_buffer_list_chisq_H1[buffer_i], search_buffer_list_chisq_L1[buffer_i], 4, 0, 0.01)


		if triggers_L1 is None:
			triggers_L1 = 0
			indiv_chisq_L1 = [1000, 0]
		if triggers_H1 is None:
			triggers_H1 = 0
			indiv_chisq_H1 = [1000,0]


		if triggers_L1 >= triggers_H1:
			triggers_out, idx_out, indiv_out, chisq_out = triggers_L1, trigger_idx_L1, indiv_L1, indiv_chisq_L1
			pivotal_ifo = 'L1'
		elif triggers_H1 > triggers_L1:
			triggers_out, idx_out, indiv_out, chisq_out= triggers_H1, trigger_idx_H1, indiv_H1, indiv_chisq_H1
			pivotal_ifo = 'H1'
		else:
			pass



		# check you haven't already counted this trigger
		if idx_out in idx_saved:
			pass

		else:
			triggers_saved.append(triggers_out)
			idx_saved.append(idx_out)
			indiv_saved.append(indiv_out)
			chisq_saved.append(chisq_out)
			cmb_chisq_saved.append(np.sqrt(0.5*(chisq_out[0]+ chisq_out[1])))

			# scrape out any particularly interesting triggers:
			if triggers_out >= loudest_thresh:
				idx_over_thresh.append(idx_out)
				chunk_over_thresh.append(chunk_i)
				buffer_no_over_thresh.append(buffer_i)
				indiv_over_thresh.append(indiv_out)
				chisq_over_thresh.append(chisq_out)


	# every n chunks, drop a new histogram file and flush the buffers
	if (chunk_i+1)%int(nchunks) == 0 or chunk_i == (min(len(chunks_L1),len(chunks_H1)) - 1):

		# could always refactor this into a nice ligolw xml

		# create a new histogram file:
		t_flag = min(chunks_L1[chunk_i].sample_times)

		# also need to add the livetime to this
		new_fname = 'zerolags_{}.txt'.format(t_flag)
		new_fname_chisq = 'chisq_{}.txt'.format(t_flag)

		# check that the background file exists
		if not os.path.exists(background_dirname):
			os.mkdir(background_dirname)

		# need to take the log10 of all the triggers
		triggers_saved_log10 = [np.log10(this_trigger) for this_trigger in triggers_saved]
		chisq_saved_log10 = [np.log10(this_chisq) for this_chisq in cmb_chisq_saved]

		hist, bin_edges = np.histogram(triggers_saved_log10, bins=hist_bins)

		chisq_hist, bin_edges_chisq = np.histogram(chisq_saved_log10, bins=hist_bins_chisq)

		if not os.path.exists(background_dirname+'bin_edges.txt'):
			np.savetxt(background_dirname+'bin_edges.txt', bin_edges)
		if not os.path.exists(background_dirname+'bin_edges_chisq.txt'):
			np.savetxt(background_dirname+'bin_edges_chisq.txt', bin_edges_chisq)

		np.savetxt(background_dirname+new_fname, hist)
		np.savetxt(background_dirname+new_fname_chisq, chisq_hist)

		print('dropped new file {}'.format(background_dirname+new_fname))

		# flush the buffer so you don't double count
		triggers_saved = []
		idx_saved = []
		idiv_saved = []
		chisq_saved = []
		cmb_chisq_saved = []

# at the very end, save the 'loudest' events so you can study them
np.savetxt(os.path.join(CWD,loudest_thresh_fname), np.column_stack([chunk_over_thresh, buffer_no_over_thresh, idx_over_thresh, indiv_over_thresh, chisq_over_thresh]))



	
