Each MATLAB file contains data from two recordings taken on the same day - one recording with orthonasal stimulation, the other with retronasal stimulation.  The two recordings are appended together into one longer data set.  They were appended so that we could do spike sorting on both files, allowing us to track single neurons across recordings.  

File naming convention: Each file name consists of a date (date of recording, MMDDYY) and experiment number.  For example 041615004_dat.mat contains data from recording number 4 on April 16, 2015.

When you open one of these files with MATLAB, you will see two variables in your workspace: 'spikes' and 'expinfo'

'spikes' is a Nx3 matrix.  Each of the N rows corresponds to one recorded spike.  Column 1 contains the spike time in sec.  Column 2 contains an arbitrary unit ID (a unique integer ID is given to each neuron).  Column 3 contains a 1 if the spike was recorded in OB, 2 if the spike was recorded in aPC.

'expinfo' is a 20x3 matrix.  Each of the 20 rows corresponds to one olfactory stimulus.  The first 10 stimuli are orthonasal.  The last 10 stimuli are retronasal.  Column 1 contains the time of the stimulus onset in sec.  Column 2 contains the odor identity (1 = 1Hexanol, 2 = Ethyl Butyrate).  Column 3 indicates which (if any) pharmacological manipulation was done to OB (0 = no drug, 1 = bicuculine, 2 = muscimol).

Stimulus is a square pulse puff of air directed through the nasal cavity with duration 1 sec.
