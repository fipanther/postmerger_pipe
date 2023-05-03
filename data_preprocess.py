"""
Load and pre-process LIGO data
"""

# want to load and process 256s segments at any one time
# can either get the data from a saved frame file, or load
# it directly using GWpy

# this should produce an appropriate data segment for processing
# and it's PSD
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py

#whittle these down later
from pycbc.types import TimeSeries
from scipy.signal.windows import tukey
import pycbc.noise
import pycbc.psd
from pycbc.filter import highpass, resample_to_delta_t, matched_filter
from pycbc.psd import interpolate, inverse_spectrum_truncation, bandlimited_interpolate


CWD = os.getcwd()

def import_data_from_file(data_dir, fname, type = 'hdf5'):

	# this will automatically throw an exception
	# if either the data_dir or the filename are not present

	# currently only support the hdf5 type:
	if type == 'hdf5':
		datafile = h5py.File(os.path.join(data_dir, fname))
	else:
		print("File type is not recognized - this version only supports HDF5!")

	# for the next bit we need to check the versioning of the h5py module
	# at some point because it will fail in the source code if they have an older
	# version of h5py

	# some of this has been borrowed from the GWOSC data tutorials and modified
	#---------------------
	# Read in strain data
	#---------------------
	# note the .value has been deprecated post h5py version 2.9.0
	strain = datafile['strain']['Strain'][:]
	ts = datafile['strain']['Strain'].attrs['Xspacing']
	sample_rate = int(1/ts)

	#-----------------------
	# We also need the metadata
	#-----------------------
	meta = datafile['meta']

	# get a 265 second working chunk
	gps_start = meta['GPSstart'][()] 
	duration = meta['Duration'][()] 
	gps_end   = gps_start + duration

	# this will make the entire timeseries from the frame file:
	
	# need to make this a PyCBC timeseries, but won't take the gps time as the start
	strain_as_timeseries = TimeSeries(strain, delta_t=ts, epoch=None)
	#return the strain, the gps start and end times 
	return strain_as_timeseries, gps_start, gps_end, duration, ts, sample_rate

def data_precondition(strain_as_timeseries, ts,
						highpass_data=True,
						highpass_freq = 20., 
						tukey_window = True, 
						tukey_rolloff = 0.5,
						Hann=False):
	
	# by default we will assume that we do not window or whiten the data
	# as these are done inside the PyCBC functions
	# If you do not use the inbuilt Welch method to do the PSD estimation 
	# then you WILL need to window the data
	if highpass_data == True:
		# this uses the inbuilt pycbc highpass method and will by default highpass
		# the data above 20Hz. You can change the highpass frequency if desired
		# might be worth changing to highpass_fir eventually for more control
		highpassed = highpass(strain_as_timeseries, highpass_freq)
		
		# snip off the spectral leakage
		#conditioned_strain  = resample_to_delta_t(conditioned_strain, 
		#	
		#								ts).crop(2, 2)

		#28/11
		# Based on Eric T's suggestion, I remove the spectral leakage with a Tukey window
		# rather than the pipeline method of just chopping off the worst offending times 
		#tukey rolloff is given in s
		window = tukey(len(highpassed), tukey_rolloff / (max(highpassed.sample_times) - min(highpassed.sample_times)))
		conditioned_strain = window * highpassed
		# end changes

	elif highpass_data == False:
		print("No highpass filter - returning unconditioned strain!")
		conditioned_strain = strain_as_timeseries


	if Hann == True:
		print("No independent Hann filter implemented yet! \n CAUTION: returning unconditioned strain!")
		conditioned_strain = strain_as_timeseries

	return conditioned_strain

def batch_strain(strain, gps_start, ts, chunk_length = 128):
	# this should return a list (for now) of chunk_length (default 128) seconds of data 
	# that can be looped through. In principle these can all be filtered completely independently
	# and then recombined at the trigger selection time
	n_chunks = int(len(strain) * ts / (chunk_length))

	# create the n chunks. Each should be overlapped eventually
	# to accomodate the windowing?

	chunks = []
	for i in range(0, n_chunks):
		this_chunk_idx_start = int(i * (chunk_length/ts))
		this_chunk_idx_end = int((i + 1) * (chunk_length/ts))
		chunk_i = strain[this_chunk_idx_start:this_chunk_idx_end]
		chunks.append(chunk_i)

	return chunks

