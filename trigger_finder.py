import numpy as np
from pycbc.types import TimeSeries

def create_search_buffers(snr_series, chisq_series, buffer_length, buffer_overlap):

	# this will create 1 second buffers in which to search for the most significant trigger

	sample_rate = int(1/snr_series.delta_t)
	
	return_buffers_snr = []
	return_buffers_chisq = []
 
	# can probably reverse engineer the buffer start/ends
	snr_length = len(snr_series)

	t_last = max(snr_series.sample_times)
	#print('tlast:{}'.format(t_last))

	# track the start index of the buffer and increment it
	this_buffer_number = 0

	# calculate the buffer length in samples:
	buffer_length_samples = int(buffer_length * sample_rate)


	next_buffer_start = int(buffer_length_samples - (buffer_overlap*sample_rate))


		# make some buffers
	while this_buffer_number <= snr_length:

		# the max length of this buffer, to stop overflows
		this_buffer_max_len = int(this_buffer_number + buffer_length_samples)
		#print('last index going into this buffer: {}'.format(this_buffer_max_len))

		# get the current second and zero-pad it if needed
		# this same logic can be used to zero-pad gaps
		if this_buffer_max_len < (snr_length):# - next_buffer_start):
			#print('in if loop, max len: {}, condition: {}'.format(this_buffer_max_len, (snr_length)))
			 #if this is true the there is enough SNR that it doesn't need zero-padding
			# ew there will be an off by one here.... 
			#print('if maximum available buffer index: {}'.format(len(snr_series)))
			this_snr_buffer = snr_series[this_buffer_number:this_buffer_max_len]
			this_chisq_buffer = chisq_series[this_buffer_number:this_buffer_max_len]
			#print('if this buffer length: {}'.format(len(this_snr_buffer)))
		# ignore this for now and come back to it
		#elif this_buffer_max_len >= (snr_length):# - next_buffer_start)
			# print('in else loop, max len: {}, condition: {}'.format(this_buffer_max_len, (snr_length)))
			# # this will be true if there is not enough SNR left to fill the whole buffer, so zero-pad
			# # always assuming here that there will be a pad at the end, in 
			# # practice might need to pad the front
			# print('else maximum available buffer index: {}'.format(len(snr_series)))
			# this_snr_buffer = snr_series[this_buffer_number:this_buffer_max_len]
			# print(min(this_snr_buffer.sample_times))
			# print('else this buffer length: {}'.format(len(this_snr_buffer)))
			# if len(this_snr_buffer) < buffer_length_samples:
			# 	# fix this
			# 	print('appending zero pad to last buffer...')
			# 	print(sample_rate - len(this_snr_buffer))
			# 	# this is currently giving a 'None' out, figure out why...
			# 	this_snr_buffer = this_snr_buffer.append_zeros(buffer_length_samples - len(this_snr_buffer))
			# 	print(this_snr_buffer)
			#print(this_snr_buffer)
			# print(max(this_snr_buffer.sample_times))
		# you made your buffer, so now increment the counter
		this_buffer_number += next_buffer_start

		return_buffers_snr.append(this_snr_buffer)
		return_buffers_chisq.append(this_chisq_buffer)



	return return_buffers_snr, return_buffers_chisq, t_last

def peaks_over_threshold(this_buffer, sngl_snr_thresh):
	# list peaks that have SNRs greater than a given SNR threshold
	peaklist = [np.abs(j)>=sngl_snr_thresh for j in this_buffer]

	return peaklist

def find_coincident_peaks(this_buffer, that_buffer, this_chisq_buffer, that_chisq_buffer, sngl_snr_thresh, cmb_snr_thresh, time_delay_max):
	"""
	Find the most significant peaks in this_buffer and find their coincident friends in 
	that_buffer
	"""

	save_peaks_this_pivotal = []
	save_peak_this_idxs = []
	save_individual = []
	save_individual_chisq = []

	# most of this will be worked in indexes rather than times
	sample_rate = int(1/this_buffer.delta_t)
	samples_in_window = int(time_delay_max*sample_rate)

	sample_start_this_buffer = int(min(this_buffer.sample_times)*sample_rate)

	this_buffer_peaklist = peaks_over_threshold(this_buffer, sngl_snr_thresh)
	
	# this gives the indicies of the True items in the peaklist
	peak_locs_this_buffer = [k for k, elem in enumerate(this_buffer_peaklist) if elem==True]

	for this_peak in peak_locs_this_buffer:
        #we can't find coincident peaks for the full window in the first or last 10ms, so chop these off
		if this_peak - samples_in_window > 0 and this_peak + samples_in_window < len(this_buffer):
			# get all the SNR timeseries from the other buffer in the coincidence window
			peaks_in_window = that_buffer[this_peak - samples_in_window:this_peak + samples_in_window]

			# pick the loudest one
			best_peak = max(np.abs(peaks_in_window))

			# NEED TO TRACK THE TIME STAMP OR INDEX HERE
			best_peak_idx = np.argmax(np.abs(peaks_in_window))

			# calculate the combined SNR of this pair
			# this is |rho|, not rho!
			combined_snr = np.sqrt(np.abs(this_buffer[this_peak])**2 + best_peak**2)

			save_peaks_this_pivotal.append(combined_snr)
			save_individual.append([np.abs(this_buffer[this_peak]), best_peak])
			save_individual_chisq.append([this_chisq_buffer[this_peak], that_chisq_buffer[best_peak_idx]])
			save_peak_this_idxs.append([this_peak, best_peak_idx+(this_peak - samples_in_window)])

	# down-select for the best of the selected peaks 
	if len(save_peaks_this_pivotal)>0:
		out = max(save_peaks_this_pivotal)
		out_idxs = save_peak_this_idxs[np.argmax(save_peaks_this_pivotal)]
		out_idxs = [idx + sample_start_this_buffer for idx in out_idxs]
		out_indiv = save_individual[np.argmax(save_peaks_this_pivotal)]
		out_indiv_chisq = save_individual_chisq[np.argmax(save_peaks_this_pivotal)]
	else:
		out = None
		out_idxs = None
		out_indiv = None
		out_indiv_chisq = None

	return out, out_idxs, out_indiv, out_indiv_chisq