clear, clc

temperature = load("./data/vostok_interpolated_temperature.csv");
age = load("./data/vostok_interpolated_ages.csv");

%% FFT

dt = 81; % sampling time, see script `./data/data_preprocessing.m`
Fs = 1/dt; % sampling frequency
L = length(temperature); % length of the signal
t = (0:L-1)*dt; % time vector
f = Fs*(0:(L/2))/L; % frequency vector

Y = fft(temperature);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

plot(1./f,P1)

%% RP und tau-RR

% recurrence threshold
epsilon = 0.08;

% MCDTS embedding gained from whole dataset analysis in `compute_embedding_parameters.jl`
embedding_taus = [0 140 46 84 112 66 57];
Yx = genembed(temperature, embedding_taus, ones(1,length(embedding_taus)));

% compute RP and tau-rr
RPx = rp(Yx,epsilon,'fan');

% tau-RR
tauRR = tau_recurrence_rate(RPx);

% Plot RP

imagesc(RPx), colormap([1 1 1; 0 0 0]), axis xy square

%% Downsampling of the tauRR
% tauRR(1:50)=0;

age_tauRR = age(1:length(tauRR));

dt_new = 4*dt; % increase sampling
new_age_tauRR = age_tauRR(1):dt_new:age_tauRR(end); % new time vector
new_tauRR = interp1(age_tauRR, tauRR, new_age_tauRR, 'pchip');

subplot(211)
plot(tauRR)
title("Original \tau-RR")
grid on
subplot(212)
plot(new_age_tauRR, new_tauRR)
title("Downsampled \tau-RR")
grid on

%% Use only the peaks of the tauRR

[pks,locs] = findpeaks(new_tauRR, 'MinPeakDistance', 30); % min peak separation of 10,000 yrs 
new_tauRR_peaks = zeros(size(new_tauRR));
new_tauRR_peaks(locs) = pks;

subplot(311)
plot(age_tauRR,tauRR)
title("original \tau-RR")
grid on
subplot(312)
plot(new_age_tauRR, new_tauRR)
title("downsampled \tau-RR")
grid on
subplot(313)
plot(new_age_tauRR, new_tauRR_peaks)
title("peaks of downsampled \tau-RR")
grid on


%% ISS of tau-RR

threshold = 0.99;
tic
[spectrum1, rho1] = inter_spike_spectrum(new_tauRR, 'threshold', threshold);
toc
tic
[spectrum2, rho2] = inter_spike_spectrum(new_tauRR_peaks, 'threshold', threshold);
toc

%% FFT of \tau-RR

Fs2 = 1/dt_new; % sampling frequency
L2 = length(new_tauRR); % length of the signal
t2 = (0:L2-1)*dt_new; % time vector
f2 = Fs2*(0:(L2/2))/L2; % frequency vector

Y2 = fft(new_tauRR);
P22 = abs(Y2/L2);
P12 = P22(1:L2/2+1);
P12(2:end-1) = 2*P12(2:end-1);

%%
subplot(411)
plot(1./f,P1)
xlim([0 200000])
title("Standard FFT of Vostok temperature")
grid on

subplot(412)
plot(new_age_tauRR(1:length(spectrum1)),spectrum1)
xlim([0 200000])
title("ISS of Vostok temperature's \tau-RR")
grid on

subplot(413)
plot(2./f2,P12)
xlim([0 200000])
title("FFT of Vostok temperature's \tau-RR")
grid on

subplot(414)
plot(new_age_tauRR(1:length(spectrum2)),spectrum2)
xlim([0 200000])
title("ISS of Vostok temperature's \tau-RR's peaks")
grid on

%%

%% Downsampling of temperature time series

dt_new = 4*dt; % increase sampling
new_age = age(1):dt_new:age(end); % new time vector
new_temperature = interp1(age, temperature, new_age, 'pchip');

subplot(211)
plot(age, temperature)
title("Original temperature")
xlim([0 400000])
grid on
subplot(212)
plot(new_age, new_temperature)
xlim([0 400000])
title("Downsampled temperature")
grid on

%% Use only the peaks of temperature time series
new_temperature = new_temperature - min(new_temperature);

[pks,locs] = findpeaks(new_temperature, 'MinPeakDistance', 30); % min peak separation of 10,000 yrs 
new_temperature_peaks = zeros(size(new_temperature));
new_temperature_peaks(locs) = pks;

subplot(311)
plot(age, temperature)
xlim([0 400000])
title("Original temperature")
grid on
subplot(312)
plot(new_age, new_temperature)
xlim([0 400000])
title("Downsampled temperature")
grid on
subplot(313)
plot(new_age, new_temperature_peaks)
xlim([0 400000])
title("Peaks of downsampled temperature")
grid on


%% ISS of temperature

threshold = 0.9;
tic
[spectrum1, rho1] = inter_spike_spectrum(new_temperature, 'threshold', threshold);
toc
tic
[spectrum2, rho2] = inter_spike_spectrum(new_temperature_peaks, 'threshold', threshold);
toc

%%
subplot(311)
plot(1./f,P1)
xlim([0 300000])
title("Standard FFT of Vostok temperature")
grid on

subplot(312)
plot(new_age_tauRR(1:length(spectrum1)),spectrum1)
xlim([0 300000])
title("ISS of downsampled Vostok temperature")
grid on

subplot(313)
plot(new_age_tauRR(1:length(spectrum2)),spectrum2)
xlim([0 300000])
title("ISS of downsampled Vostok temperature's peaks")
grid on

%%

save("./data/new_tauRR.csv", "new_tauRR", "-ascii")
save("./data/new_age_tauRR.csv", "new_age_tauRR", "-ascii")