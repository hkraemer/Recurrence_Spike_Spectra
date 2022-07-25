% In this script we 1) downsample the original data and 2) compute
% FFT-spectrograms on these downsampled time series.

clear, clc

f_ce = load("./data/frequencies_ce.mat").frequencies_ce;
f_gb = load("./data/frequencies_gb.mat").frequencies_gb;

%% Downsample to a sampling time of one minute/ 20 sec / 30 sec

T_ce = 0.2; % original sampling time in secs
T_gb = 1;

T_ce_ds = 20; % new sampling time 20 s
T_gb_ds = 20;

factor_ce = 1/T_ce * T_ce_ds;
factor_gb = 1/T_gb * T_gb_ds;

f_ce_ds = f_ce(1:factor_ce:end);
f_gb_ds = f_gb(1:factor_gb:end);

L_ce_ds = length(f_ce_ds);
L_gb_ds = length(f_gb_ds);

t_ce = (0:L_ce_ds-1)*T_ce_ds; % time vector
t_gb = (0:L_gb_ds-1)*T_gb_ds;


%% Short time FFT

T_ce = T_ce_ds; % sampling time in secs
T_gb = T_gb_ds;
Fs_ce = 1/T_ce; % sampling frequency
Fs_gb = 1/T_gb;

L_ce = length(f_ce); % length of signal
L_gb = length(f_gb);

L_block_ce = 4320; % block length covering one full day (24 hours)
L_block_gb = 4320;

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

%% Save results 

% save("./data/frequencies_ce_downsampled_20s.csv", "f_ce_ds", "-ascii")
% save("./data/frequencies_gb_downsampled_20s.csv", "f_gb_ds", "-ascii")

save("./data/fft_spec_ce.csv","spec_ce","-ascii")
save("./data/fft_spec_gb.csv","spec_gb","-ascii")


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

