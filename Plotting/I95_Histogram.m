% Histogram representing 
% the periods of the first 20 pulses when dynamic entrainment is applied.
% I0=95 and g_Ca=4 is used.

clear all
close all
clc

% Gap data
gaps = [24.42, 28.03, 25.43, 26.7, 26, 26.04, 26.03, 26.06, 26.02, 26.04, ...
        26.03, 26.06, 26.01, 26.07, 26.02, 26.06, 26.04, 26.05, 26.02, 26.04]; % Amp=12, L=65

% gaps = [20.29, 21.49, 19.74, 20.38, 19.91, 20.1, 19.92, 20.09, 19.93, 20.92, ...
%     19.97, 20, 19.94, 20.01, 19.94, 19.99, 19.95, 20.01, 19.98, 19.97]; % Amp=12, L=70

% gaps = [51.35,55.58,55.58,55.58,55.58,55.58,55.58,55.58,55.58,55.58,55.58, ...
%     55.58,55.58,55.58,55.58,55.58,55.58,55.58,55.58,55.58]; % Amp=12, L=35


% Create histogram
figure(1)
% histogram(gaps, 'BinWidth', 0.5); 
histogram(gaps, 'BinWidth', 0.5, 'Normalization', 'pdf');
xlabel('Gap between successive pulses',Interpreter='latex');
ylabel('Density',Interpreter='latex');
grid on;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',18);

% period data
periods=[89.42, 93.03, 90.43, 91.7, 91, 91.04, 91.03, 91.06, 91.02, 91.04, ...
         91.03, 91.06, 91.01, 91.07, 91.02, 91.06, 91.04, 91.05, 91.02, 91.04];% Amp=12, L=65

% periods=[90.29, 91.49, 89.74, 90.38, 89.91, 90.1, 89.92, 90.09, 89.93, 90.92, ...
%     89.97, 90, 89.94, 90.01, 89.94, 89.99, 89.95, 90.01, 89.98, 89.97];% Amp=12, L=70 

% periods=[86.35,90.58,90.58,90.58,90.58,90.58,90.58,90.58,90.58,90.58, ...
%     90.58,90.58,90.58,90.58,90.58,90.58,90.58,90.58,90.58,90.58];% Amp=12, L=35 

% periods=[97.67,145.6,147.35,144.41,144.36,144.36,144.36,144.36, ...
%     144.36,144.36,144.36,144.36,144.36,144.36,144.36,144.36,...
%     144.36,144.36,144.36,144.36];% Amp=12, L=95 - didn't give 1:1 

figure(2)
histogram(periods, 'BinWidth', 0.5, 'Normalization', 'pdf');
xlabel('$T$',Interpreter='latex');
ylabel('Density',Interpreter='latex');
grid on;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',18);