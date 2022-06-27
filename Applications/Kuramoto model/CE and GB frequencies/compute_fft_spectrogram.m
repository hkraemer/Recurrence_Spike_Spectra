% FFT of the frequency-data

clear, clc

f_ce = load("./data/frequencies_ce.mat").frequencies_ce;
f_gb = load("./data/frequencies_gb.mat").frequencies_gb;

%%

L_ce = length(f_ce); % length of the signal
L_gb = length(f_gb);
T_ce = 0.2; % sampling time in secs
T_gb = 1;
Fs_ce = 1/T_ce; % sampling frequency
Fs_gb = 1/T_gb;

t_ce = (0:L_ce-1)*T_ce; % time vector
t_gb = (0:L_gb-1)*T_gb;

% FFT
Y_ce = fft(f_ce(1:L_ce));
Y_gb = fft(f_gb(1:L_gb));

P2_ce = abs(Y_ce/L_ce);
P1_ce = P2_ce(1:L_ce/2+1);
P1_ce(2:end-1) = 2*P1_ce(2:end-1);
fs_ce= Fs_ce*(0:(L_ce/2))/L_ce;

P2_gb = abs(Y_gb/L_gb);
P1_gb = P2_gb(1:L_gb/2+1);
P1_gb(2:end-1) = 2*P1_gb(2:end-1);
fs_gb= Fs_gb*(0:(L_gb/2))/L_gb;

%%

fs_ce= (1/(0.2/60))*(0:(L_ce/2))/L_ce;
fs_gb= (1/(1/60))*(0:(L_gb/2))/L_gb;


figure
subplot(211)
semilogx(fs_ce,P1_ce)
title('Single-Sided Amplitude Spectrum of F_{CE}')
xlabel('f (Hz)')
ylabel('|P1(f)|')
xlim([10e-5 1])
grid on

subplot(212)
semilogx(fs_gb,P1_gb)
title('Single-Sided Amplitude Spectrum of F_{GB}')
xlabel('f (Hz)')
ylabel('|P1(f)|')
xlim([10e-5 1])
grid on



