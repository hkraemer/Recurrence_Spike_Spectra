% For the CENOGRID dataset we perform a windowed analysis with the spike 
% powerspectrum, in order to obtain a wavelet-like spectrum. Here we
% analyze the eccentricity time series from Laskar 2011.
 
% Here: TRAJECTORY BASED ON MCDTS embedding, FNN + continuity statistic (script 
% `compute_embedding_parameters.jl` in this folder)

clear, clc

data = load("../data/Laskar/ecc_Laskar.txt");

% sampling time is dt = 5,000 yrs

%% Set parameters

epsilon = 0.1; % recurrence threshold
windowsize = 1200; % windowsize for trajectory
ws = 1; % windowstep
N = length(data); % length of the time series
rho_thres = 0.9; % thresholds for inter spike spectrum
M = length(1:ws:N-windowsize); % number of computed spectra

%% Compute inter spike spectrogram in windowed analysis

window_tau = 200; % sampling time is 5.000 yrs so this corresponds to 1 Mio yrs
length_of_spectra = ceil(window_tau/2);

spectrum_ecc_mcdts = zeros(length_of_spectra,M); 
spectrum_ecc_ne = zeros(length_of_spectra,M);

taus_mcdts = [0 8 62 115 4 3]; % from script `compute_embedding_parameters.jl`

%%
cnt = 1;
tic
for i = 1:ws:M
    display(i)
    % extract time series
    x = data(i:i+windowsize);
    % 1) MCDTS embedding gained from whole dataset analysis
%     Yx = genembed(x,taus_mcdts, ones(1,length(taus_mcdts)));

    % compute RP and tau-rr
%     RPx = rp(Yx,epsilon,'var');
    RPx_ne = rp(x,epsilon,'var');
%     tauRRx = tau_recurrence_rate(RPx);
    tauRRx_ne = tau_recurrence_rate(RPx_ne);
    % compute spike powerspectra
%     spectrum_ecc_mcdts(:,cnt) = inter_spike_spectrum(tauRRx(1:window_tau), "method", "STLS", "threshold", rho_thres);
    spectrum_ecc_ne(:,cnt) = inter_spike_spectrum(tauRRx_ne(1:window_tau), "method", "STLS", "threshold", rho_thres);

    cnt = cnt + 1;
    
end
toc

% save("./computed spectra/spectrum_ecc_mcdts.mat", "spectrum_ecc_mcdts")
save("./computed spectra/spectrum_ecc_ne.mat", "spectrum_ecc_ne")
