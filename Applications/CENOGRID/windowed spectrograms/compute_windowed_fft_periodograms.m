% FFT - ANALYSIS
%
% For the CENOGRID dataset we perform a windowed FFT analysis (short time 
% Fourier Transform with a Hamming window), in order to obtain a wavelet-like 
% spectrum.

clear, clc

data = load("../data/detrended.txt");

t = data(:,1);
O18_raw = flipud(data(:,2));
C13_raw = flipud(data(:,3));
t = flipud(t);
t = t.*1000000;

O18 = load("../data/O18_filtered_time_reverse.txt");
C13 = load("../data/C13_filtered_time_reverse.txt");

ecc = load("../data/Laskar/ecc_Laskar.txt");


%% Compute windowed periodogramm

dt = abs(mean(diff(t))); % Sampling period
Fs = 1/dt; % Sampling frequency

%%
windowsize2 = 1000; % windowsize for trajectory
windowsize1 = 100; % windowsize for the spectrum
ws = 1; % windowstep


[s, f, ~] = spectrogram(O18_raw,windowsize2,windowsize2-ws,[],1);
pxx_raw = abs(s).^2;
[s, ~, ~] = spectrogram(O18,windowsize2,windowsize2-ws,[],1);
pxx_filtered = abs(s).^2;
[s, ~, ~] = spectrogram(C13_raw,windowsize2,windowsize2-ws,[],1);
pxy_raw = abs(s).^2;
[s, ~, ~] = spectrogram(C13,windowsize2,windowsize2-ws,[],1);
pxy_filtered = abs(s).^2;
[s, ~, ~] = spectrogram(ecc,windowsize2,windowsize2-ws,[],1);
p_ecc = abs(s).^2;

f = f .* Fs;


save("periodograms_raw.mat",'f','pxx_raw','pxy_raw','p_ecc')
save("periodograms_filtered.mat",'f','pxx_filtered','pxy_filtered','p_ecc')
