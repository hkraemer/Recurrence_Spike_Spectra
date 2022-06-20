% read data

clear; close all; clc; 
addpath('raw_data');

filename = 'raw_data/monofrequent_sound.wav'; % simplest sound
filename = 'raw_data/many_harmonics.wav';  % nonlinear vibration with higher harmonics
filename = 'raw_data/multi_freq_sounds.wav'; % most complex sound

[y, fs] = audioread(filename);
t = [0:length(y)-1] ./ fs;

figure; 
plot(t, y); xlabel('time'); ylabel('sound pressure level [Pa]'); 