# Signal_processing - Signal recovery techniques (linear distortion)
The project involves passing an audio recording through several typical signal processing stages, including signal sampling, filtering, digital filter implementation, and investigating signals in both the time and frequency domains. Audacity was used to generate the audio file.

The Matlab code has been structured into tasks as you will notice in the source file and below are some other important comments:

# Task C
The resampling factor is fe'/fe = 160/441. Since this is not an integer value, two stages are needed for this subtask: one for interpolation to increase the sampling frequency and one for decimation to decrease the same frequency. The order of these operations is not chosen randomly. Thus, interpolation is limited to the old value of fe/2, and decimation is limited to the new value of the sampling frequency. A band limitation is needed to avoid spectral overlap, which is accomplished using an IIR filter.

# Task D
The role of the IIR filter is to prevent relevant spectral overlap in the frequency domain. Thus, this filter is of high order for strong attenuation (difficult to build), requiring signal band limitation up to fe/2. Regarding the 7 kHz cutoff frequency, this value is the upper limit of the human voice frequency range (vocal band ϵ [50 Hz and 7 kHz]). This 7 kHz value also ensures compliance with Nyquist's condition (the maximum frequency of the signal could be 8 kHz = fe/2).
