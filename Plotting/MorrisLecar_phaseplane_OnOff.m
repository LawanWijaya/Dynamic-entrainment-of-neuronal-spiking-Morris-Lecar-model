% on and off pulse
% red and black
% nullclines

clear all;
close all;
clc
% These values are from Ermentrout Table 3.1 first column (Hopf)
Amp=25; % No square pulse because Amp=0
php=0.5;
L=30;
per=2*L;
% I0= 95; % Change Accordingly
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

%time
start_pulse=0;
end_pulse=3;
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

I=I0+ Amp*(heaviside(mod(t,per)-per*(1-php))); % square pulse

% Initialize arrays for different forcing states
V_forcing_on = V; V_forcing_on(I==I0) = nan;  % Red segment
V_forcing_off = V; V_forcing_off(I>I0) = nan; % Black segment
n_forcing_on = n; n_forcing_on(I==I0) = nan;
n_forcing_off = n; n_forcing_off(I>I0) = nan;


%% Plotting
% % % figure(1);
% % % subplot(4,1,[1 3])
% % % % plot(t-tstart,ySol(:,1), '-',Color=[0 0 1])
% % % plot(t-tstart, V_forcing_off, 'k', 'LineWidth', 2.5); % Plot V in black when I is 0
% % % hold on
% % % plot(t-tstart, V_forcing_on, 'r', 'LineWidth', 2.5); % Plot V in red when I is non-zero
% % % hold on
% % % % ylabel('Membrane Voltage',Interpreter='latex')
% % % ylabel('$V$',Interpreter='latex')
% % % % xlabel('$t$',Interpreter='latex')
% % % set(groot,'defaultAxesTickLabelInterpreter','latex');
% % % set(groot,'defaulttextinterpreter','latex');
% % % set(groot,'defaultLegendInterpreter','latex');
% % % set(gca,'TickLabelInterpreter','latex','fontsize',16);
% % % xticklabels([]); % Remove x-axis labels in this subplot
% % % xlim([0 tend-tstart]);
% % % 
% % % subplot(4,1,4)
% % % plot(t - tstart, I, '-r', 'LineWidth', 2.5)
% % % xlabel('$t$', 'Interpreter', 'latex')
% % % ylabel('$Amp$', 'Interpreter', 'latex')
% % % set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
% % % xlim([0 tend - tstart]);
% % % hold on;
% % % 
% % % % Detect where forcing is OFF
% % % tolerance = 1e-3;
% % % zeroSegments = abs(I - I0) < tolerance;
% % % 
% % % % Find segment boundaries
% % % edges = diff([0; zeroSegments; 0]);
% % % startIdx = find(edges == 1);
% % % endIdx = find(edges == -1) - 1;
% % % 
% % % % Plot black overlays for forcing off
% % % for k = 1:length(startIdx)
% % %     idx_range = startIdx(k):endIdx(k);
% % %     plot(t(idx_range) - tstart, I(idx_range), 'k', 'LineWidth', 2);
% % % end

figure(2)
hold on;
% fimplicit(dVdt,[-70 40 -0.1 0.5],Color='r'); % V-nullcline
fimplicit(dVdt_0,[-70 50 -0.1 0.6],'LineWidth',2,'Color','k','LineStyle','--'); % V-nullcline when forcing off
fimplicit(dVdt_A,[-70 50 -0.1 0.6],'LineWidth',2,'Color','r','LineStyle','--'); % V-nullcline when forcing on
xlim([-70 50])
ylim([-0.1 0.6])
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
box


function dYdt = Morris(t,y, Amp, php, per, Iapp, phi, gCa, V1, V2, V3, V4, ECa, EK, EL, gK, gL, Cm)

% I=Iapp+ Amp*(1-heaviside(mod(t,per)-per(1-php)));
I=Iapp+ Amp*(heaviside(mod(t,per)-per*(1-php)));

V=y(1);
n=y(2);

mV=0.5*(1+tanh((V-V1)/V2));
tV=1/cosh((V-V3)/(2*V4));
n1V=0.5*(1+tanh((V-V3)/V4));

dVdt=(I-gL*(V-EL)-gK*n*(V-EK)-gCa*mV*(V-ECa))/Cm;
dndt=phi*(n1V-n)/tV;
dYdt=[dVdt;dndt;];
end