% We prepare the temperature Vostok data (interpolation)

clear, clc

temperature = load("vostok_temperature.csv");
age = load("vostok_age.csv");

% due to the sampling differences focus on the first 400000 yrs
idx = find(age>400000);
age = age(1:idx(1));
temperature = temperature(1:idx(1));

% mean sampling time
dt = floor(mean(abs(diff(age))));

% make new time vector according to the mean sampling time
new_age = age(1):dt:age(end);

% interpolate
new_temperature = interp1(age, temperature, new_age,'pchip');
new_temperature = new_temperature(2:end);
new_age = new_age(2:end);

%%

subplot(211)
plot(age,temperature,'.-')
grid on
subplot(212)
plot(new_age,new_temperature,'.-')
grid on

%% Save results
save("vostok_interpolated_temperature.csv", "new_temperature", "-ascii")
save("vostok_interpolated_ages.csv", "new_age", "-ascii")
