% We plot the downsampled time series and its autocorrelation 
% (computed in `compute_fft_spectrogram_downsampled`), the FFT-spectrograms 
% (computed in `compute_fft_spectrogram_downsampled`) and the according ISS 
% (computed in `compute_iss_spectrogram_downsampled.jl`).

clear, clc

%% Load the computed data

% (downsampled) time series (computed in `compute_fft_spectrogram_downsampled`)
f_ce = load("./data/frequencies_ce_downsampled_20s.csv");
f_gb = load("./data/frequencies_gb_downsampled_20s.csv");

% FFT spec (computed in `compute_fft_spectrogram_downsampled`)
spec_ce = load("./data/fft_spec_ce.csv");
spec_gb = load("./data/fft_spec_gb.csv");

% ISS spec (computed in `compute_iss_spectrogram_downsampled.jl`)
iss_spec_ce = load("./data/iss_spec_ce_standard_20_0_95_0_8.csv");
iss_spec_gb = load("./data/iss_spec_gb_standard_20_0_95_0_8.csv");

%% Set time and frequency vectors

dt_s = 20; % sampling time is sec
dt_m = 1/ (60/dt_s);

% time/period vector for ISS 
N = length(iss_spec_ce);
t = (1:N)*dt_m;

% time/period/frequency vector for FFT spectrogram
Fs_ce = 60/dt_s; % sampling frequency min^-1
Fs_gb = 60/dt_s;

L_block_ce = 4320; % block length covering one full day (24 hours, 60/dt_s * 60 * 24))
L_block_gb = 4320;

t_ce = (0:L_block_ce-1)*dt_m; % time vector for block in min
t_gb = (0:L_block_gb-1)*dt_m;

fs_ce= Fs_ce*(0:(L_block_ce/2))/L_block_ce; % frequency vector in min^1
fs_gb= Fs_gb*(0:(L_block_gb/2))/L_block_gb;

%% Plotting

font = 22; % fontsize
lw1 = 2.5;
lw2 = 2;

figure('Units','normalized','Position',[.001 .001 .999 .999])
subplot(221)
semilogx(1./fs_ce, spec_ce, 'LineWidth', lw1)
title('Single-Sided Amplitude Spectrum of F_{CE}')
xline(15,'r--','LineWidth', lw1)
xline(30,'r--','LineWidth', lw1)
xline(60,'r--','LineWidth', lw1)
% xlabel('period (min)')
ylabel('|power/period|')
xlim([1 100])
set(gca, 'LineWidth', lw2, 'FontSize', font)
grid on


subplot(222)
semilogx(1./fs_gb, spec_gb, 'LineWidth', lw1)
title('Single-Sided Amplitude Spectrum of F_{GB}')
xline(15,'r--','LineWidth', lw1)
xline(30,'r--','LineWidth', lw1)
xline(60,'r--','LineWidth', lw1)
ylabel('|power/period|')
xlim([1 100])
set(gca, 'LineWidth', lw2, 'FontSize', font)
grid on

subplot(223)
semilogx(t, iss_spec_ce(1:300), 'LineWidth', lw1)
title('ISS of \tau-RR of F_{CE}')
xline(15,'r--','LineWidth', lw1)
xline(30,'r--','LineWidth', lw1)
xline(60,'r--','LineWidth', lw1)
ylabel('|power/period|')
xlabel('period (min)')
xlim([1 100])
ylim([0 0.01])
set(gca, 'LineWidth', lw2, 'FontSize', font)
grid on

subplot(224)
semilogx(t, iss_spec_gb(1:300), 'LineWidth', lw1)
title('ISS of \tau-RR of F_{GB}')
xline(15,'r--','LineWidth', lw1)
xline(30,'r--','LineWidth', lw1)
xline(60,'r--','LineWidth', lw1)
ylabel('|power/period|')
xlabel('period (min)')
xlim([1 100])
ylim([0 0.01])
set(gca, 'LineWidth', lw2, 'FontSize', font)
grid on

%% Plot time series and autocorrelation

font = 24; % fontsize
lw1 = 3;
lw2 = 2;

init = 4459; 
factor = 3; % corresponds to 24/factor hours
t_len = L_block_ce/factor + 1; 

% compute autocorrelation
t_len2 = 300;
ac_ce = acorr(f_ce, t_len2, 0);
ac_gb = acorr(f_gb, t_len2, 0);

%%

figure('Units','normalized','Position',[.001 .001 .999 .999])
subplot(211)
plot(t_ce(1:t_len), f_ce(init: init+L_block_ce/factor), '.-', 'LineWidth', lw1)
title('F_{CE}')
ylabel('frequency (Hz)')
xlabel('time (min)')
xlim([0 24/factor*60])
set(gca, 'LineWidth', lw2, 'FontSize', font)
grid on

subplot(212)
plot(t_ce(1:t_len2+1), ac_ce(1:end,2), '.-', 'LineWidth', lw1)
ylabel('normal. acorr')
xlabel('time lag \tau (min)')
xlim([1 100])
set(gca, 'LineWidth', lw2, 'FontSize', font)
grid on

figure('Units','normalized','Position',[.001 .001 .999 .999])
subplot(211)
plot(t_gb(1:t_len), f_gb(init: init+L_block_gb/factor),  '.-', 'LineWidth', lw1)
title('F_{GB}')
xlabel('time (min)')
ylabel('frequency (Hz)')
xlim([0 24/factor*60])
set(gca, 'LineWidth', lw2, 'FontSize', font)
grid on

subplot(212)
plot(t_ce(1:t_len2+1), ac_gb(1:end,2), '.-', 'LineWidth', lw1)
ylabel('normal. acorr')
xlabel('time lag \tau (min)')
xlim([1 100])
set(gca, 'LineWidth', lw2, 'FontSize', font)
grid on
