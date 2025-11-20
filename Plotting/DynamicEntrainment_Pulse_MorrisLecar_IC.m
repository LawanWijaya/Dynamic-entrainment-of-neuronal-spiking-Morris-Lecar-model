% Using the pulse predicted by using dynamic entrainment
% Initial conditions are determined by running MorrisLecar_for_IC.m

% using initial_values from Morris_Lecar_forIC.m
clear all;
close all;
clc
% These values are from Ermentrout Table 3.1 first column (Hopf)

% From generated dataset:
% Lowest Amp = 10.00103421,  highesr amp = 79.99996823
% Lowest L = 9.96494535,  highest L = 99.99973864

Amp=15; % No square pulse because Amp=0
L=15;
per=2*L;
php=0.5;
I0= 80; % Change Accordingly
phi=0.04;
gCa=4.4;
V1=-1.2;
V2=18;
V3=2;
V4=30;
ECa=120;
EK=-84;
EL=-60;
gK=8;
gL=2;
Cm=20;

% pulse_starts array represent the pulse onsets determined by dynamic entrainment

% pulse_starts = [20,103.71, 187.47, 271.09, 354.83, 438.45];%Amp=30, L=20 % starting time of each pulse
% pulse_starts = [20,98.71, 177.28, 256.1, 334.78, 413, 492, 571, 650,729,808,887,966]; % starting time of each pulse
% pulse_starts=[20,98,176,254,332,410,488,566,644];
% pulse_starts=[20, 171.18,322,473,624,775];
% pulse_starts=[20, 108.5, 197, 285.5, 374, 462.5, 551, 639, 727];
% pulse_starts=[20, 108.43, 197.11, 285.49, 373.74, 462.37]; % Amp=25, L=35
% pulse_starts=[20, 122.83,226,329,432,535];
pulse_starts=[20, 122.83, 224.95, 327.82, 429.95, 532.82]; % Amp=15, L=15
% pulse_starts=[20, 107.49 195.29 282.64 370.33 457.67 ]; % Amp=25, L=30
% pulse_starts = [20,109.42,212.25,315.08,417.91,520.74,623.57]; % Amp=15,L=15, t0 and start is tnext



tstart=0;
tend=pulse_starts(end)+L+40;
tIni=0;
time_range=tIni:0.01:tend;
time_range2=tstart:0.01:tend;



%% ODE Solver
% % % initial_values=[0,0.1];
options=odeset('AbsTol',1.e-8,'relTol',1.e-9,'InitialStep',1.e-3,'MaxStep',0.01);
% % % [t1, ySol1]=ode15s(@Morris, time_range, initial_values, options, ...
% % %     Amp, I0, phi, gCa, V1, V2, V3, V4, ECa, EK, EL, gK, gL, Cm, ...
% % %     L, pulse_starts);

% % % initial_values=[ySol1(end,1),ySol1(end,2)];
% initial_values=[-24.3681,0.1207]; %Amp=15, L=15
% initial_values=[8.9288,0.1751]; %Amp=25, L=30
% initial_values=[36.3847,0.2680]; %Amp=30, L=20
% initial_values=[23.8606,0.4093]; %Amp=25, L=30, pulse on from t0=20
initial_values=[ -26.3832,0.1212]; %Amp=15, L=15, pulse on from t0=20


[t, ySol]=ode15s(@Morris, time_range2, initial_values, options, ...
    Amp, I0, phi, gCa, V1, V2, V3, V4, ECa, EK, EL, gK, gL, Cm, ...
    L, pulse_starts);

% Solution of the differential equations
V=ySol(:,1); %y(1)
n=ySol(:,2); %y(2)

%% For V-nullcline 
% When forcing off
F=0;
dVdt_0 = @(V,n) I0+F-gL*(V-EL)-gK*n*(V-EK)-gCa*(0.5*(1+tanh((V-V1)/V2)))*(V-ECa);

% When forcing on
F=Amp;
dVdt_A = @(V,n) I0+F-gL*(V-EL)-gK*n*(V-EK)-gCa*(0.5*(1+tanh((V-V1)/V2)))*(V-ECa);

%% For n-nullcline
% Not affected by forcing
dndt=@(V,n) 0.5*(1+tanh((V-V3)/V4))-n;

%% Pulse informed by machine learning data
I = I0 * ones(size(t));


for i = 1:length(pulse_starts)
    pulse_on = t >= pulse_starts(i) & t < (pulse_starts(i) + L);
    I(pulse_on) = I(pulse_on) + Amp;
end



%% Plotting
figure(1);
subplot(4,1,[1 3])
plot(t-tstart, V, 'b', 'LineWidth', 2.5); % Plot V in black when I is 0
hold on
ylabel('$V$',Interpreter='latex')
% xlabel('$t$',Interpreter='latex')
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',18);
xticklabels([]); % Remove x-axis labels in this subplot
xlim([0 tend-tstart]);

subplot(4,1,4)
plot(t - tstart, I, '-r', 'LineWidth', 2.5)
xlabel('$t$', 'Interpreter', 'latex')
ylabel('$I_{app}$', 'Interpreter', 'latex')
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
xlim([0 tend - tstart]);
hold on;


figure(2)
hold on;
% fimplicit(dVdt,[-70 40 -0.1 0.5],Color='r'); % V-nullcline
fimplicit(dVdt_0,[-90 50 -0.2 1],'LineWidth',2,'Color','k','LineStyle','--'); % V-nullcline when forcing off
fimplicit(dVdt_A,[-90 50 -0.2 1],'LineWidth',2,'Color','r','LineStyle','--'); % V-nullcline when forcing on

hold on;
fimplicit(dndt,'LineWidth',2,'Color','g'); %n-nullcline
plot(V, n, 'b', 'LineWidth', 2.5); % Plot Forcing On in Red
xlabel('$V$',Interpreter='latex')
ylabel('$n$',Interpreter='latex')
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',16);
xlim([-80 50])
ylim([-0.2 0.6])
box

function dYdt = Morris(t,y, Amp, I0, phi, gCa, V1, V2, V3, V4, ECa, EK, EL, gK, gL, Cm, ...
    L, pulse_starts)


%% Pulse informed by machine learning data
I = I0 * ones(size(t));


for i = 1:length(pulse_starts)
    pulse_on = t >= pulse_starts(i) & t < (pulse_starts(i) + L);
    I(pulse_on) = I(pulse_on) + Amp;
end
%%
V=y(1);
n=y(2);

mV=0.5*(1+tanh((V-V1)/V2));
tV=1/cosh((V-V3)/(2*V4));
n1V=0.5*(1+tanh((V-V3)/V4));

dVdt=(I-gL*(V-EL)-gK*n*(V-EK)-gCa*mV*(V-ECa))/Cm;
dndt=phi*(n1V-n)/tV;
dYdt=[dVdt;dndt;];
end