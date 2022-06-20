% time integrations from the initial conditions on all solutions

clear; close all; clc; 


%% Reference damping

% Parameter
m = 1;
mu = 0.65; 
kx = 11.2; ky = 20; kz1 = 100; kz2 = 20; 
klin1 = 10; klin2 = 10; knl1 = 5; knl2 = 5; 
cx = 0.02; cy = 0.02; cz1 = 0.02; cz2 = 0.02; clin1 = 0.02; clin2 = 0.02; 

% Parametervektor für die mitgelieferte ODE-Funktion
p = [ m, mu, kx, ky, kz1, kz2, klin1, klin2, knl1, knl2, cx, cy, cz1, cz2, clin1, clin2];

% wichtig!!  kx = 11.2;

y01 = [10.9743095964043;-17.2723918339629;-0.715936487677040;-5.12073221742318;-3.92807376346864;7.47143811701081;-2.16006060628278;0.580684576350545];
y02 = [14.8920132488256;-21.8396202917156;-3.10740915202540;-6.89290790235976;7.43511993609540;-6.20257134522162;-6.70046783700049;-5.35710902370515];
y03 = [18.3237595308644;-22.2963061108126;-3.12114019866061;-4.42314721558563;8.40978666208839;-3.21915875635150;-11.3637582770204;-3.86338826491144];
y04 = [21.3896856643288;-23.9787884941944;-2.77469389852530;-3.09661478233620;6.11696381130031;7.90622513101437;-20.4529708534388;-3.17556005203984];
y05 =[18.7525762111720;-19.8871056755027;-3.25757215892014;-2.46100000573605;-57.6690791349349;65.3639305196251;1.47570036661630;5.33242536681206];

% time integration which should arrive at a LC solution
tspan = [0, 800]; 
odefun = matcont_model_4DOF;

p(3) = 11.2; % just to make sure
par = num2cell(p);

% do the time integration
[T1,Y1] = ode45(@(t,y) odefun{2}(t,y,par{:}),tspan,y01);
[ pos_max, ~, ~, ~] = computeEnvelopeTimeSeries(Y1(:,1));
y_period1 = Y1(pos_max(end-1):pos_max(end),:);

[T2,Y2] = ode45(@(t,y) odefun{2}(t,y,par{:}),tspan,y02);
[ pos_max, ~, ~, ~] = computeEnvelopeTimeSeries(Y2(:,1));
y_period2 = Y2(pos_max(end-1):pos_max(end),:);

[T3,Y3] = ode45(@(t,y) odefun{2}(t,y,par{:}),tspan,y03);
[ pos_max, ~, ~, ~] = computeEnvelopeTimeSeries(Y3(:,1));
y_period3 = Y3(pos_max(end-1):pos_max(end),:);

[T4,Y4] = ode45(@(t,y) odefun{2}(t,y,par{:}),tspan,y04);
[ pos_max, ~, ~, ~] = computeEnvelopeTimeSeries(Y4(:,1));
y_period4 = Y4(pos_max(end-1):pos_max(end),:);

[T5,Y5] = ode45(@(t,y) odefun{2}(t,y,par{:}),tspan,y05);
[ pos_max, ~, ~, ~] = computeEnvelopeTimeSeries(Y5(:,1));
y_period5 = Y5(pos_max(end-1):pos_max(end),:);


figure; 
ax1 = subplot(5,1,1);plot(T1, Y1(:,1)); hold on; 
ylabel('$x$', 'interpreter', 'latex'); set(gca, 'Xtick',[]); legend('starting from LC1');

ax2 = subplot(5,1,2);plot(T2, Y2(:,1));
ylabel('$x$', 'interpreter', 'latex'); set(gca, 'Xtick',[]); legend('starting from LC2');

ax3 = subplot(5,1,3);plot(T3, Y3(:,1));
ylabel('$x$', 'interpreter', 'latex'); set(gca, 'Xtick',[]); legend('starting from LC3');

ax4 = subplot(5,1,4);plot(T4, Y4(:,1));
ylabel('$x$', 'interpreter', 'latex'); set(gca, 'Xtick',[]); legend('starting from LC4');

ax5 = subplot(5,1,5);plot(T5, Y5(:,1));
ylabel('$x$', 'interpreter', 'latex'); xlabel('time'); legend('starting from LC5');
legend('starting from LC1', 'starting from LC2', 'starting from LC3', 'starting from LC4', 'starting from LC5')
linkaxes([ax1, ax2, ax3, ax4, ax5],'xy'); 


figure; 
plot(y_period1(:,1), y_period1(:,5), '-^'); hold on;
plot(y_period2(:,1), y_period2(:,5), '-');
plot(y_period3(:,1), y_period3(:,5), '--s');
plot(y_period4(:,1), y_period4(:,5), '-x');
plot(y_period5(:,1), y_period5(:,5), '--o');hold on;
xlabel('$x$', 'interpreter', 'latex'); ylabel('$\dot{x}$', 'interpreter', 'latex');
l = legend('starting from LC1', 'starting from LC2', 'starting from LC3', 'starting from LC4', 'starting from LC5'); 
l.Location='northeastoutside';


%% Weakly damped configuration


clear; close all; clc; 

% Parameter
m = 1;
mu = 0.65; 
kx = 11; ky = 20; kz1 = 100; kz2 = 20; 
klin1 = 10; klin2 = 10; knl1 = 5; knl2 = 5; 
cx = 0.002; cy = 0.002; cz1 = 0.002; cz2 = 0.002; clin1 = 0.002; clin2 = 0.002; 

% Parametervektor für die mitgelieferte ODE-Funktion
p = [ m, mu, kx, ky, kz1, kz2, klin1, klin2, knl1, knl2, cx, cy, cz1, cz2, clin1, clin2];

% wichtig!!  kx = 11;

y01 = [11.1721211017737;-17.3690860702920;-0.686972834157252;-5.03769796694696;-0.362275788544908;0.584973464686077;-0.141188656739585;0.0250350526099915];
y02 = [15.2891024338951;-21.9755054959772;-3.12117904938034;-6.73453493776826;-1.71899933822516;1.61183524377635;0.760639249264874;0.502403950808016];
y03 = [18.4179541501575;-22.4492402183431;-3.25704504463500;-4.58132315348899;-1.97456727287208;1.84907348864090;0.662746648913814;0.322805931920670];
y04 = [12.9264985617386;-12.6266064120359;-2.85972884783172;-1.36968632435559;-92.6743236084569;101.317433261554;7.15703313142275;9.01821845045663]; 
y05 = [21.8536169036225;-24.1123243407510;-3.38740376529686;-3.25381616486731;-2.51800585401895;2.84109055475120;0.0113365772194647;0.198817016377782]; 
 

% time integration which should arrive at a LC solution
tspan = [0:0.00001:500]; 
odefun = matcont_model_4DOF;

p(3) = 11; % just to make sure
par = num2cell(p);

% do the time integration
[T1,Y1] = ode45(@(t,y) odefun{2}(t,y,par{:}),tspan,y01);
[ pos_max1, ~, ~, ~] = computeEnvelopeTimeSeries(Y1(:,1));
y_period1 = Y1(pos_max1(end-1):pos_max1(end),:);

[T2,Y2] = ode45(@(t,y) odefun{2}(t,y,par{:}),tspan,y02);
[ pos_max2, ~, ~, ~] = computeEnvelopeTimeSeries(Y2(:,1));
y_period2 = Y2(pos_max2(end-1):pos_max2(end),:);

[T3,Y3] = ode45(@(t,y) odefun{2}(t,y,par{:}),tspan,y03);
[ pos_max3, ~, ~, ~] = computeEnvelopeTimeSeries(Y3(:,1));
y_period3 = Y3(pos_max3(end-1):pos_max3(end),:);

[T4,Y4] = ode45(@(t,y) odefun{2}(t,y,par{:}),tspan,y04);
[ pos_max4, ~, ~, ~] = computeEnvelopeTimeSeries(Y4(:,1));
y_period4 = Y4(pos_max4(end-1):pos_max4(end),:);

[T5,Y5] = ode45(@(t,y) odefun{2}(t,y,par{:}),tspan,y05);
[ pos_max5, ~, ~, ~] = computeEnvelopeTimeSeries(Y5(:,1));
y_period5 = Y5(pos_max5(end-1):pos_max5(end),:);


figure; 
step = 10;
ax1 = subplot(5,1,1);plot(T1(1:step:end), Y1(1:step:end,1)); hold on; 
ylabel('$x$', 'interpreter', 'latex'); set(gca, 'Xtick',[]); legend('starting from LC1');

ax2 = subplot(5,1,2);plot(T2(1:step:end), Y2(1:step:end,1));
ylabel('$x$', 'interpreter', 'latex'); set(gca, 'Xtick',[]); legend('starting from LC2');

ax3 = subplot(5,1,3);plot(T3(1:step:end), Y3(1:step:end,1));
ylabel('$x$', 'interpreter', 'latex'); set(gca, 'Xtick',[]); legend('starting from LC3');

ax4 = subplot(5,1,4);plot(T4(1:step:end), Y4(1:step:end,1));
ylabel('$x$', 'interpreter', 'latex'); set(gca, 'Xtick',[]); legend('starting from LC4');

ax5 = subplot(5,1,5);plot(T5(1:step:end), Y5(1:step:end,1));
ylabel('$x$', 'interpreter', 'latex'); xlabel('time'); legend('starting from LC5');
xlim([0, 200]); 
linkaxes([ax1, ax2, ax3, ax4, ax5],'xy');

pos_max1(Y1(pos_max1,1)<0)=[];
pos_max2(Y2(pos_max2,1)<0)=[];
pos_max3(Y3(pos_max3,1)<0)=[];
pos_max4(Y4(pos_max4,1)<0)=[];
pos_max5(Y5(pos_max5,1)<0)=[];
figure; 
plot(T1(pos_max1), Y1(pos_max1,1), 'LineWidth',0.8); hold on; 
plot(T2(pos_max2), Y2(pos_max2,1), 'LineWidth',0.8);
plot(T3(pos_max3), Y3(pos_max3,1), 'LineWidth',0.8);
plot(T4(pos_max4), Y4(pos_max4,1), 'LineWidth',0.8);
plot(T5(pos_max5), Y5(pos_max5,1), 'LineWidth',0.8);
xlabel('time'); ylabel('upper envelope $x$', 'interpreter', 'latex');

figure; 
plot(y_period1(:,1), y_period1(:,5), '-^'); hold on;
plot(y_period2(:,1), y_period2(:,5), '-');
plot(y_period3(:,1), y_period3(:,5), '--s');
plot(y_period4(:,1), y_period4(:,5), '-x');
plot(y_period5(:,1), y_period5(:,5), '--o');hold on;
xlabel('$x$', 'interpreter', 'latex'); ylabel('$\dot{x}$', 'interpreter', 'latex');
l = legend('LC1', 'LC2', 'LC3', 'LC4', 'LC5');
l.Location='northeastoutside'; 

