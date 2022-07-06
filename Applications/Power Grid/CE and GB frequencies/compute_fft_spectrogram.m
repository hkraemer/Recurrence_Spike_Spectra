% FFT of the frequency-data according to the paper Meyer, Anvari & Kantz
% 2020, Chaos(30)

clear, clc

f_ce = load("./data/frequencies_ce.mat").frequencies_ce;
f_gb = load("./data/frequencies_gb.mat").frequencies_gb;

%% Short time FFT

T_ce = 0.2; % sampling time in secs
T_gb = 1;
Fs_ce = 1/T_ce; % sampling frequency
Fs_gb = 1/T_gb;

L_ce = length(f_ce); % length of signal
L_gb = length(f_gb);

L_block_ce = 432000; % block length covering one full day (24 hours)
L_block_gb = 86400;

t_ce = (0:L_block_ce-1)*T_ce; % time vector
t_gb = (0:L_block_gb-1)*T_gb;

fs_ce= (1/(T_ce))*(0:(L_block_ce/2))*60/L_block_ce; % frequency vector
fs_gb= (1/(T_gb))*(0:(L_block_gb/2))*60/L_block_gb;

tic
spectrum_ce = fft_spec(f_ce, L_block_ce);
toc
tic
spectrum_gb = fft_spec(f_gb, L_block_gb);
toc

%% Averaging over all windows

spec_ce = mean(spectrum_ce,2);
spec_gb = mean(spectrum_gb,2);

%% Plotting

figure
subplot(211)
semilogx(1./fs_ce, spec_ce)
title('Single-Sided Amplitude Spectrum of F_{CE}')
xlabel('period (min)')
ylabel('|P1(f)|')
grid on


subplot(212)
semilogx(1./fs_gb, spec_gb)
title('Single-Sided Amplitude Spectrum of F_{GB}')
xlabel('period (min)')
ylabel('|P1(f)|')
grid on


%%
function spectrum = fft_spec(x, blocksize)
    T = length(x);
    blocks = floor(T/blocksize);
    
    spectrum = zeros(blocksize/2+1, blocks);
    for i = 1:blocks
        time = 1+(i-1)*blocksize:(i*blocksize);
        Y = fft(x(time));
        P2 = abs(Y/blocksize);
        P1 = P2(1:blocksize/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        spectrum(:,i) = P1;
    end
end

