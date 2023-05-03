This repository contains code that was used in Panther & Lasky 2023 (arXiv:2303.10847)

The main program is main_fd.py - this will generate a simple postmerger template waveform and perform (FD) matched filtering between GW detector strain and the relevant template. It then performs a two-detector coincidence search for 'significant' (signal-like) triggers. Note that all the candidates it identifies are zerolags - i.e. technically if you want to do a real search and estimate significance, you should perform timeslides to identify the distribution of background triggers associated with each zerolag. The main program draws functions from the pre-processing package (data_preprocess.py) and the trigger finder (trigger_finder.py). 
Details of the algorithm can be found in the above paper, or email fiona.panther@uwa.edu.au
The output of main_fd.py is a series of text files that contain binned combined SNRs and combined chisq for each zerolag, as well as a list of the most significant triggers. The histogram generator script allows the visualizations of the zerolag histograms. 

TODO: add script that allows one to parse and investigte the most 'significant' triggers output by the main_fd.py script  
