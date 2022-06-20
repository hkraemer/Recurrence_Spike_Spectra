% read data

clear; close all; clc;
addpath('raw_data');

filename = 'raw_data/Signal_test_n_1.mat';
% filename = 'raw_data/Signal_test_n_2.mat';
% filename = 'raw_data/Signal_test_n_3.mat';
% filename = 'raw_data/Signal_test_n_4.mat';
% filename = 'raw_data/Signal_test_n_5.mat';

load(filename);
y = Signal_test_n_1; clear Signal_test_n_1

t = y(:,1); y(:,1) = [];
fs = 1/(t(2)-t(1)); % 100kHz sampling rate


figure;
leg_entries = {'fixed block velocity', 'moving block acceleration', 'tangential force', 'normal force'};
ax = cell(4,1);
for i=1:4
    ax{i,1} = subplot(4,1,i);
    plot(t, y(:,i));
    legend(leg_entries{i});
end
xlabel('time');
linkaxes([ax{:}], 'x');