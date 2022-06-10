% For the CENOGRID dataset we perform a windowed analysis with the spike 
% powerspectrum, in order to obtain a wavelet-like spectrum. Here we
% analyze the CENOGRID time series.
 
% Here: TRAJECTORY BASED ON MCDTS embedding, FNN + continuity statistic (script 
% `compute_embedding_parameters.jl` in this folder)

clear, clc

data = load("../data/detrended.txt");

t = data(:,1);
t = flipud(t);
t = t.*1000000;

O18 = load("../data/O18_filtered_time_reverse.txt");
C13 = load("../data/C13_filtered_time_reverse.txt");

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

spectrum_O18_mcdts = zeros(length_of_spectra,M); 
spectrum_O18_ne = zeros(length_of_spectra,M);
spectrum_C13_mcdts = zeros(length_of_spectra,M); 
spectrum_C13_ne = zeros(length_of_spectra,M);

taus_O18_mcdts = [0  50  104  35  43  149  39  47  45]; % from script `compute_embedding_parameters.jl`
taus_C13_mcdts = [0  136  32  123  130  126  133  128]; % from script `compute_embedding_parameters.jl`

%%
cnt = 1;
tic
for i = 1:ws:M
    display(i)
    % extract time series
    x = O18(i:i+windowsize);
    y = C13(i:i+windowsize);
    % 1) MCDTS embedding gained from whole dataset analysis
    Yx = genembed(x,taus_O18_mcdts, ones(1,length(taus_O18_mcdts)));
    Yy = genembed(x,taus_C13_mcdts, ones(1,length(taus_C13_mcdts)));

    % compute RP and tau-rr
    RPx = rp(Yx,epsilon,'var');
    RPx_ne = rp(x,epsilon,'var');
    RPy = rp(Yy,epsilon,'var');
    RPy_ne = rp(y,epsilon,'var');
    tauRRx = tau_recurrence_rate(RPx);
    tauRRx_ne = tau_recurrence_rate(RPx_ne);
    tauRRy = tau_recurrence_rate(RPy);
    tauRRy_ne = tau_recurrence_rate(RPy_ne);
    % compute spike powerspectra 
    spectrum_O18_mcdts(:,cnt) = inter_spike_spectrum(tauRRx(1:window_tau), "method", "STLS", "threshold", rho_thres);
    spectrum_O18_ne(:,cnt) = inter_spike_spectrum(tauRRx_ne(1:window_tau), "method", "STLS", "threshold", rho_thres);
    spectrum_C13_mcdts(:,cnt) = inter_spike_spectrum(tauRRy(1:window_tau), "method", "STLS", "threshold", rho_thres);
    spectrum_C13_ne(:,cnt) = inter_spike_spectrum(tauRRy_ne(1:window_tau), "method", "STLS", "threshold", rho_thres);
    
    cnt = cnt + 1;
    
end
toc

save("./computed spectra/spectrum_O18_mcdts.mat", "spectrum_O18_mcdts")
save("./computed spectra/spectrum_O18_ne.mat", "spectrum_O18_ne")
save("./computed spectra/spectrum_C13_mcdts.mat", "spectrum_C13_mcdts")
save("./computed spectra/spectrum_C13_ne.mat", "spectrum_C13_ne")
