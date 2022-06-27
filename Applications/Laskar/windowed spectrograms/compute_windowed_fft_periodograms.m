% FFT - ANALYSIS
%
% For the LASKAR dataset we perform a windowed FFT analysis (short time 
% Fourier Transform with a Hamming window), in order to obtain a wavelet-like 
% spectrum of the eccentricity time series.

clear, clc

ecc = load("../data/Laskar/ecc_Laskar.txt");

dt = 5000; % Sampling period
Fs = 1/dt; % Sampling frequency

t = linspace(0,6.7100e+07,length(ecc)); % time vector


%% Compute windowed periodogramm

windowsize2 = 1000; % windowsize for trajectory
ws = 1; % windowstep

[s, f, ~] = spectrogram(ecc,windowsize2,windowsize2-ws,[],1);
p_ecc = abs(s).^2;

f = f .* Fs;

save("./computed spectra/periodograms_short_time_fourier.mat", 'f', 'p_ecc')
