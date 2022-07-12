clear, clc

f_ce = load("./data/frequencies_ce.mat").frequencies_ce;
f_gb = load("./data/frequencies_gb.mat").frequencies_gb;

%% Downsample to a sampling time of one minute

T_ce = 0.2; % sampling time in secs
T_gb = 1;

T_ce_ds = 20; % new sampling time
T_gb_ds = 20;

factor_ce = 1/T_ce * T_ce_ds;
factor_gb = 1/T_gb * T_gb_ds;

f_ce_ds = f_ce(1:factor_ce:end);
f_gb_ds = f_gb(1:factor_gb:end);

L_ce_ds = length(f_ce_ds);
L_gb_ds = length(f_gb_ds);

t_ce = (0:L_ce_ds-1)*T_ce_ds; % time vector
t_gb = (0:L_gb_ds-1)*T_gb_ds;

% % checking
% subplot(211)
% plot(t_ce(1:21000), f_ce(1:21000), '.-')
% grid on
% subplot(212)
% plot(t_ce_ds(1:21000/factor_ce), f_ce_ds(1:21000/factor_ce), '.-')
% grid on

%% Short time FFT

T_ce = T_ce_ds; % sampling time in secs
T_gb = T_gb_ds;
Fs_ce = 1/T_ce; % sampling frequency
Fs_gb = 1/T_gb;

L_ce = length(f_ce); % length of signal
L_gb = length(f_gb);

L_block_ce = 1440; % block length covering one full day (24 hours)
L_block_gb = 1440;

t_ce = (0:L_block_ce-1)*T_ce; % time vector for block
t_gb = (0:L_block_gb-1)*T_gb;

fs_ce= Fs_ce*(0:(L_block_ce/2))*60/L_block_ce; % frequency vector
fs_gb= Fs_gb*(0:(L_block_gb/2))*60/L_block_gb;

tic
spectrum_ce = fft_spec(f_ce_ds, L_block_ce);
toc
tic
spectrum_gb = fft_spec(f_gb_ds, L_block_gb);
toc

%% Averaging over all windows

spec_ce = mean(spectrum_ce,2);
spec_gb = mean(spectrum_gb,2);

%% Plotting

% iss_spec_ce = load("./data/iss_spec_ce_cao_0_8_0_5.csv");
% iss_spec_gb = load("./data/iss_spec_gb_cao_0_8_0_5.csv");
% iss_spec_ce = load("./data/iss_spec_ce_hegger_0_8_0_5.csv");
% iss_spec_gb = load("./data/iss_spec_gb_hegger_0_8_0_5.csv");
% iss_spec_ce = load("./data/iss_spec_ce_hegger_0_9_0_5.csv");
% iss_spec_gb = load("./data/iss_spec_gb_hegger_0_9_0_5.csv");
% iss_spec_ce = load("./data/iss_spec_ce_standard_60_0_9_0_5.csv");
% iss_spec_gb = load("./data/iss_spec_gb_standard_60_0_9_0_5.csv");
% iss_spec_ce = load("./data/iss_spec_ce_cao_60_0_95_0_8.csv");
% iss_spec_gb = load("./data/iss_spec_gb_cao_60_0_95_0_8.csv");
% iss_spec_ce = load("./data/iss_spec_ce_cao_60_0_99_0_8.csv");
% iss_spec_gb = load("./data/iss_spec_gb_cao_60_0_99_0_8.csv");
% iss_spec_ce = load("./data/iss_spec_ce_cao_STLS_60_0_95_0_8.csv");
% iss_spec_gb = load("./data/iss_spec_gb_cao_STLS_60_0_95_0_8.csv");
% iss_spec_ce = load("./data/iss_spec_ce_standard_60_0_95_0_8.csv");
% iss_spec_gb = load("./data/iss_spec_gb_standard_60_0_95_0_8.csv");
iss_spec_ce = load("./data/iss_spec_ce_standard_20_0_95_0_8.csv");
iss_spec_gb = load("./data/iss_spec_gb_standard_20_0_95_0_8.csv");

%%
t = (1:300)*0.333333333333333;

figure
subplot(221)
semilogx(1./fs_ce, spec_ce)
title('Single-Sided Amplitude Spectrum of F_{CE}')
xline(15,'r--')
xline(30,'r--')
xline(60,'r--')
xlabel('period (min)')
ylabel('|P1(f)|')
xlim([1 100])
grid on


subplot(222)
semilogx(1./fs_gb, spec_gb)
title('Single-Sided Amplitude Spectrum of F_{GB}')
xline(15,'r--')
xline(30,'r--')
xline(60,'r--')
xlabel('period (min)')
ylabel('|P1(f)|')
xlim([1 100])
grid on

subplot(223)
semilogx(t, iss_spec_ce(1:300))
title('ISS of F_{CE}')
xline(15,'r--')
xline(30,'r--')
xline(60,'r--')
xlabel('period (min)')
ylabel('|P1(f)|')
xlim([1 100])
ylim([0 0.01])
grid on

subplot(224)
semilogx(t, iss_spec_gb(1:300))
title('ISS of F_{GB}')
xline(15,'r--')
xline(30,'r--')
xline(60,'r--')
xlabel('period (min)')
ylabel('|P1(f)|')
xlim([1 100])
ylim([0 0.01])
grid on

%% ISS for downsampled time series

x = f_ce_ds;
blocksize = L_block_ce;

T = length(x);
blocks = floor(T/blocksize);

spectrum = zeros(blocksize/2+1, blocks);

i = 1;

time = 1+(i-1)*blocksize:(i*blocksize);
data = x(time);

[Y, TAU_VALS, TS_VALS] = pecuzal_embedding(data);

%%
save("frequencies_ce_downsampled_20s.csv", "f_ce_ds", "-ascii")
save("frequencies_gb_downsampled_20s.csv", "f_gb_ds", "-ascii")

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

