% generate data

clear; close all; clc; 
addpath('computations');

m = 1;
mu = 0.65; 
kx = 11; ky = 20; kz1 = 100; kz2 = 20; 
klin1 = 10; klin2 = 10; knl1 = 5; knl2 = 5; 
cx = 0.002; cy = 0.002; cz1 = 0.002; cz2 = 0.002; clin1 = 0.002; clin2 = 0.002; 

% Parametervektor für die mitgelieferte ODE-Funktion
p = [ m, mu, kx, ky, kz1, kz2, klin1, klin2, knl1, knl2, cx, cy, cz1, cz2, clin1, clin2];

% initial conditions for the different co-existing solutions (no 4 is the
% chaotic one)
% y01 = [11.1721211017737;-17.3690860702920;-0.686972834157252;-5.03769796694696;-0.362275788544908;0.584973464686077;-0.141188656739585;0.0250350526099915];
% y02 = [15.2891024338951;-21.9755054959772;-3.12117904938034;-6.73453493776826;-1.71899933822516;1.61183524377635;0.760639249264874;0.502403950808016];
% y03 = [18.4179541501575;-22.4492402183431;-3.25704504463500;-4.58132315348899;-1.97456727287208;1.84907348864090;0.662746648913814;0.322805931920670];
y04 = [12.9264985617386;-12.6266064120359;-2.85972884783172;-1.36968632435559;-92.6743236084569;101.317433261554;7.15703313142275;9.01821845045663]; 
% y05 = [21.8536169036225;-24.1123243407510;-3.38740376529686;-3.25381616486731;-2.51800585401895;2.84109055475120;0.0113365772194647;0.198817016377782]; 
 
% time integration which should arrive at a LC solution
fs = 1000;
tspan = [0:1/fs:500]; 
odefun = matcont_model_4DOF;

p(3) = 11; % just to make sure that we have the correct parameter set here
par = num2cell(p);

% do the time integration
disp('running time integration ...');
[t,y] = ode45(@(t,y) odefun{2}(t,y,par{:}),tspan,y04);
disp('finished');

figure; 
n = size(y,2); % number of states
ax = cell(n,1);
for i=1:n
    ax{i,1} = subplot(n, 1, i);
    plot(t, y(:,i)); 
    ylabel(['state ', num2str(i)]);
end
xlabel('time'); 
linkaxes([ax{:}], 'x');


