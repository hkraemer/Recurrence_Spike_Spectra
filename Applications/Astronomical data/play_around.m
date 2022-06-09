clear, clc

%% Load the data

data1 = load("./data/0791780601_imaging_S001");
data2 = load("./data/0670920501_timing_S003"); % this here takes ages for the spike spectrum and also yields rand deficient warnings

% data1 = load("./data/0411083201_image_S600.txt");
% data2 = load("./data/0502030101_timing_S003.txt");

% data1 = load("./data/0153950701_imaging_S005.txt");
% data2 = load("./data/0153951201_timing_S005.txt");

subplot(211)
plot(data1(:,1),data1(:,2))
grid on
subplot(212)
plot(data2(:,1),data2(:,2))
grid on

%% Construct RP

threshold = 0.10;

RP = rp(data2(:,2),threshold,'var');
tauRR = tau_recurrence_rate(RP);
spectrum = inter_spike_spectrum(tauRR, 'method', "STLS");

figure
subplot(4,1,[1,2])
imagesc(RP), colormap([1 1 1; 0 0 0]), axis xy square
subplot(4,1,3)
plot(tauRR)
grid on
subplot(4,1,4)
plot(spectrum)
grid on

%%

save("test_tauRR.csv", "tauRR", "-ascii")