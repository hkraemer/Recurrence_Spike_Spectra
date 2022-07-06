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

L_block = 512; % block length

t_ce = (0:L_block-1)*T_ce; % time vector
t_gb = (0:L_block-1)*T_gb;

fs_ce= (1/(T_ce/60))*(0:(L_block/2))/L_block; % frequency vector
fs_gb= (1/(T_gb/60))*(0:(L_block/2))/L_block;

blocks_ce = floor(L_ce/L_block);
blocks_gb = floor(L_gb/L_block);

% tic
% spectrum_ce = fft_spec(f_ce, L_block);
% toc
% tic
% spectrum_gb = fft_spec(f_gb, L_block);
% toc

tic
[spectrum_ce, ff_ce, ts_ce] = stft(f_ce, Fs_ce, 'OverlapLength', 0, 'FFTLength',512);
toc
tic
[spectrum_gb, ff_gb, ts_gb] = stft(f_gb, Fs_gb, 'OverlapLength', 0, 'FFTLength',512);
toc

%%
spec_ce = mean(spectrum_ce,2);
spec_gb = mean(spectrum_gb,2);

%% Plotting
i = 100000;

P2 = abs(spectrum_ce(:,i)/L_block);
P1 = P2(1:L_block/2+1);
P1(2:end-1) = 2*P1(2:end-1);

fs_ce = Fs_ce*(0:(L_block/2))/L_block;
fs_gb = Fs_gb*(0:(L_block/2))/L_block;

figure
subplot(211)
semilogx(fs_ce, P1 )
title('Single-Sided Amplitude Spectrum of F_{CE}')
xlabel('f (Hz)')
ylabel('|P1(f)|')
% xlim([10e-5 1])
grid on

P2 = abs(spectrum_gb(:,i)/L_block);
P1 = P2(1:L_block/2+1);
P1(2:end-1) = 2*P1(2:end-1);

subplot(212)
semilogx(fs_gb, P1)
title('Single-Sided Amplitude Spectrum of F_{GB}')
xlabel('f (Hz)')
ylabel('|P1(f)|')
% xlim([10e-5 1])
grid on

% P2 = abs(spec_ce/L_block);
% P1 = P2(1:L_block/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% 
% figure
% subplot(211)
% semilogx(fs_ce, P1)
% title('Single-Sided Amplitude Spectrum of F_{CE}')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
% % xlim([10e-5 1])
% grid on
% 
% P2 = abs(spec_gb/L_block);
% P1 = P2(1:L_block/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% 
% subplot(212)
% semilogx(fs_gb, P1)
% title('Single-Sided Amplitude Spectrum of F_{GB}')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
% % xlim([10e-5 1])
% grid on

%%


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

