% on and off pulse
% red and black
% Similar to MachineLearning_Pulse_MorrisLecar.m 
% but 
% limit cycle is red and black based on forcing being on and off.
% starts from the steady state values instead of given initial conditions

clear all;
close all;
clc

% These values are from Ermentrout Table 3.1 first column (Hopf)
% Book: Mathematical Foundations of Neuroscience by G. Bard Ermentrout and David H. Terman

% From generated dataset:
% Lowest Amp = 10.00103421,  highesr amp = 79.99996823
% Lowest L = 9.96494535,  highest L = 99.99973864

Amp=25; % No square pulse because Amp=0
L=35;
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

% starting time of each pulse
% pulse_starts = [20,103.71, 187.47, 271.09, 354.83]; % Amp=20, L=20. t0=20
% pulse_starts=[20, 108.43, 197.11, 285.49, 373.74, 462.37]; % Amp=25, L=35
% pulse_starts=[20, 108.43]; % Amp=25, L=35
% pulse_starts=[20, 122.83, 224.95, 327.82, 429.95, 532.82]; % Amp=15, L=15
% pulse_starts=[20, 122.83, 224.95]; % Amp=15, L=15
% pulse_starts=[20, 107.49 195.29 282.64 370.33 457.67 ]; % Amp=25, L=30
pulse_starts=[20, 107.49]; % Amp=25, L=30


% pulse_starts = [20,98.71, 177.28, 256.1, 334.78, 413, 492, 571, 650,729,808,887,966]; % Amp=10, L=20, t0=20
% pulse_starts=[20,98,176,254,332,410,488,566,644];
% pulse_starts=[20, 171.18,322,473,624,775];
tstart=100; % where we start
tend=pulse_starts(end)+L+60;
tIni=0;
time_range=tIni:0.01:tend;
time_range2=tstart:0.01:tend;


%% ODE Solver
initial_values=[0,0.1];
options=odeset('AbsTol',1.e-8,'relTol',1.e-9,'InitialStep',1.e-3,'MaxStep',0.01);
[t1, ySol1]=ode15s(@Morris, time_range, initial_values, options, ...
    Amp, I0, phi, gCa, V1, V2, V3, V4, ECa, EK, EL, gK, gL, Cm, ...
    L, pulse_starts);

initial_values=[ySol1(end,1),ySol1(end,2)];

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

% Initialize arrays for different forcing states
V_forcing_on = V; V_forcing_on(I==I0) = nan;  % Red segment
V_forcing_off = V; V_forcing_off(I>I0) = nan; % Black segment
n_forcing_on = n; n_forcing_on(I==I0) = nan;
n_forcing_off = n; n_forcing_off(I>I0) = nan;


%% Plotting

% Time series diagram
% % % figure(1);
% % % subplot(4,1,[1 3])
% % % plot(t-tstart, V_forcing_off, 'k', 'LineWidth', 2.5); % Plot V in black when I is 0
% % % hold on
% % % plot(t-tstart, V_forcing_on, 'r', 'LineWidth', 2.5); % Plot V in red when I is non-zero
% % % ylabel('$V$',Interpreter='latex')
% % % % xlabel('$t$',Interpreter='latex')
% % % set(groot,'defaultAxesTickLabelInterpreter','latex');
% % % set(groot,'defaulttextinterpreter','latex');
% % % set(groot,'defaultLegendInterpreter','latex');
% % % set(gca,'TickLabelInterpreter','latex','fontsize',22);
% % % xticklabels([]); % Remove x-axis labels in this subplot
% % % xlim([0 tend-tstart]);
% % % 
% % % subplot(4,1,4)
% % % plot(t - tstart, I, '-r', 'LineWidth', 2.5)
% % % xlabel('$t$', 'Interpreter', 'latex')
% % % ylabel('$I_{app}$', 'Interpreter', 'latex')
% % % set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 22);
% % % xlim([0 tend - tstart]);
% % % hold on;


figure(2)
hold on;
% fimplicit(dVdt,[-70 40 -0.1 0.5],Color='r'); % V-nullcline
fimplicit(dVdt_0,[-90 50 -0.2 1],'LineWidth',2,'Color','k','LineStyle','--'); % V-nullcline when forcing off
fimplicit(dVdt_A,[-90 50 -0.2 1],'LineWidth',2,'Color','r','LineStyle','--'); % V-nullcline when forcing on

hold on;
fimplicit(dndt,'LineWidth',2,'Color','g'); %n-nullcline
plot(V_forcing_on, n_forcing_on, 'r', 'LineWidth', 2.5); % Plot Forcing On in Red
plot(V_forcing_off, n_forcing_off, 'k', 'LineWidth', 2.5); % Plot Forcing Off in Black
xlabel('$V$',Interpreter='latex')
ylabel('$n$',Interpreter='latex')
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',22);
xlim([-70 50])
ylim([-0.1 0.6])
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