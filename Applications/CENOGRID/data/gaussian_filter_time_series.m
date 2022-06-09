% apply a Gaussian filter to the detrended time series.

clear, clc

data = load("../data/detrended.txt");

t = data(:,1);
O18 = data(:,2);
C13 = data(:,3);

t = flipud(t);
O18 = flipud(O18);
C13 = flipud(C13);

Filterlength = 10;

w = gausswin(Filterlength);
w = w/sum(w);

O18_filtered_time_reverse = filter(w,1,O18);
C13_filtered_time_reverse = filter(w,1,C13);

%% 
save("O18_filtered_time_reverse.txt", "O18_filtered_time_reverse","-ascii")
save("C13_filtered_time_reverse.txt", "C13_filtered_time_reverse","-ascii")
