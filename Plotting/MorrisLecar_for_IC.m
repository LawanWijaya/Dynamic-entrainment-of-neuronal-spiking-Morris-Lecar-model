% on and off pulse
% red and black
% nullclines

% Morris Lecar model with a pulse of 50% duty cycle
% Solve this to get initial conditions for
% DynamicEntrainment_Pulse_MorrisLecar.m

clear all;
close all;
clc
% These values are from Ermentrout Table 3.1 first column (Hopf)

% From generated dataset:
% Lowest Amp = 10.00103421,  highesr amp = 79.99996823
% Lowest L = 9.96494535,  highest L = 99.99973864

Amp=30; % No square pulse because Amp=0
L=20;
per=2*L;
php=0.5;
I0= 80; % Change Accordingly
phi=0.04;
gCa=4.4; % From Ermentrout
% gCa=4;% From Khan FHN paper
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


%time
start_pulse=0;
end_pulse=10;
tIni=0;
tstart=per*start_pulse;
tend=per*end_pulse;
time_range=tIni:0.01:tend;
time_range2=tstart:0.01:tend;

%% ODE Solver
initial_values=[0,0.1];
options=odeset('AbsTol',1.e-8,'relTol',1.e-9,'InitialStep',1.e-3,'MaxStep',0.01);
[t1, ySol1]=ode15s(@Morris, time_range, initial_values, options, ...
    Amp, php, per, I0, phi, gCa, V1, V2, V3, V4, ECa, EK, EL, gK, gL, Cm);

initial_values=[ySol1(end,1),ySol1(end,2)];

[t, ySol]=ode15s(@Morris, time_range2, initial_values, options, ...
    Amp, php, per, I0, phi, gCa, V1, V2, V3, V4, ECa, EK, EL, gK, gL, Cm);

% Solution of the differential equations
V=ySol(:,1); %y(1)
n=ySol(:,2); %y(2)

% I=I0+ Amp*(heaviside(mod(t,per)-per*(1-php))); % square pulse

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

%% Give red or black to limit cycle based on F =A or F=0

% % % I=I0+ Amp*(heaviside(mod(t,per)-per*(1-php))); % square pulse
I=I0+ Amp*(heaviside(-mod(t-20,per)+per*(1-php)).* (t>=20)); % square pulse, start from on


%% Plotting
figure(1);
subplot(4,1,[1 3])
plot(t-tstart, V, 'b', 'LineWidth', 2.5); % Plot V in black when I is 0
hold on
ylabel('$V$',Interpreter='latex')
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
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
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
set(gca,'TickLabelInterpreter','latex','fontsize',18);
xlim([-80 50])
ylim([-0.2 0.6])
box

fprintf("Initial Values for Dynamic Entrainment Pulse: [%f %f]",initial_values(1),initial_values(2));

function dYdt = Morris(t,y, Amp, php, per, I0, phi, gCa, V1, V2, V3, V4, ECa, EK, EL, gK, gL, Cm)

% I=Iapp+ Amp*(1-heaviside(mod(t,per)-per(1-php)));
% % % I=I0+ Amp*(heaviside(mod(t,per)-per*(1-php)));
I=I0+ Amp*(heaviside(-mod(t-20,per)+per*(1-php)).* (t>=20)); % square pulse, start from on

V=y(1);
n=y(2);

mV=0.5*(1+tanh((V-V1)/V2));
tV=1/cosh((V-V3)/(2*V4));
n1V=0.5*(1+tanh((V-V3)/V4));

dVdt=(I-gL*(V-EL)-gK*n*(V-EK)-gCa*mV*(V-ECa))/Cm;
dndt=phi*(n1V-n)/tV;
dYdt=[dVdt;dndt;];
end