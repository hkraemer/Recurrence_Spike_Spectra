% For the CENOGRID dataset we perform a windowed analysis with the spike 
% powerspectrum, in order to obtain a wavelet-like spectrum. Here we
% analyze the eccentricity time series from Laskar 2011.
 
% Here: TRAJECTORY BASED ON PECUZAL embedding (script 
% `compute_embedding_parameters.jl` in this folder)

clear, clc

data = load("../data/Laskar/ecc_Laskar.txt");


%% Set parameters

epsilon = 0.05; % recurrence threshold
windowsize2 = 1000; % windowsize for trajectory
windowsize1 = 100; % windowsize for the spectrum
ws = 1; % windowstep
N = length(data); % length of the time series
taus = 0:150; % delay range passed to PECUZAL
rho_thres = 0.95; % threshold for inter spike spectrum

M = length(1:ws:N-windowsize2);

%% Compute inter spike spectrogram in windowed analysis

spectrum_ecc_pec = zeros(windowsize1,M); 
spectrum_ecc_ne = zeros(windowsize1,M); 

cnt = 1;
tic
for i = 1:ws:M
    display(i)
    % extract time series
    x = data(i:i+windowsize2);
    % 1) Pecuzal embedding gained from whole dataset analysis
    Yx = genembed(x,[0, 8, 4, 21, 15, 33, 27, 128, 74, 135, 64], ones(1,10));

    % compute RP and tau-rr
    RPx = rp(Yx,epsilon,'var');
    tauRRx = tau_recurrence_rate(RPx);
    % compute spike powerspectra
    spectrum_ecc_pec(:,cnt) = get_spike_spectrum(tauRRx(1:windowsize1), lambda);

   
    cnt = cnt + 1;
    
end
toc


save("spectrum_ecc_pec.mat", "spectrum_ecc_pec")


%% Compute spike powerspectra in windowed analysis

spectrum_ecc_pec2 = zeros(windowsize1,M); 

cnt = 1;
tic
for i = 1:ws:M
    display(i)
    % extract time series
    x = data(i:i+windowsize2);
    % 1) Pecuzal embedding gained from whole dataset analysis
    Yx = genembed(x,[0,13,7,10,56,27,3,5,77,42,20,17,15], ones(1,13));

    % compute RP and tau-rr
    RPx = rp(Yx,epsilon,'var');
    tauRRx = tau_recurrence_rate(RPx);
    % compute spike powerspectra
    spectrum_ecc_pec2(:,cnt) = get_spike_spectrum(tauRRx(1:windowsize1), lambda);

   
    cnt = cnt + 1;
    
end
toc


save("spectrum_ecc_pec2.mat", "spectrum_ecc_pec2")