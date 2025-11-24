% Instead of adding a third pulse
% we are changing the initial conditions: IC = [V,n] and starting time t0

% commented item 1 from the original code 
% (because we want system start to start from IC instead of the steady state)

% commented item 3 because that is related to the second pulse (stage 3 and
% 4)

% The responses from this code is input to the machine learning model to
% calculate pulse timing.

function [TRAINING_DATASET_V1, TRAINING_DATASET_V2] = MorrisLecar_After_First_Pulse(seed, COUNTER)

% MorrisLecar_After_First_Pulse - Generate data based on the Morris–Lecar model

% using MATLAB's parallel computing capabilities.

%

% Usage:

%   [DS1, DS2] = MorrisLecar_After_First_Pulse(seed, COUNTER)

%

% Inputs:

%   seed    - Random seed for reproducibility.

%   COUNTER - Number of simulations (dataset lines) to run.

%

% Outputs:

%   TRAINING_DATASET_V1 - Array with features [Amp, L, t1, p1, t2, p2, tnextmax, t3max, p3max]

%   TRAINING_DATASET_V2 - Array with features [Amp, L, L1, p1, L2, p2, tnextmax, l3max, p3max]



% Model Parameters for Morris-Lecar Model

V1 = -1.2; V2 = 18; V3 = 2; V4 = 30;

gL = 2; EL = -60; gK = 8; EK = -84; gCa = 4.4; ECa = 120;

phi = 0.04; C = 20; php = 0.5;

stepsize = 1; lw = 3;

% Stimulation parameter bounds
% L_l = 50; L_u = 65;
% A_l = 95; A_u = 105;

A_l = 35;
L_l = 15;
I0=80;

% Use provided COUNTER or default to 1 if not given

if nargin < 2

    COUNTER = 1;

end



% Set the random seed for reproducibility.

% rng(seed, 'twister');

% Preallocate output arrays.

% % % % % TRAINING_DATASET_V1 = zeros(COUNTER, 9);
% % % % % 
% % % % % TRAINING_DATASET_V2 = zeros(COUNTER, 9);
TRAINING_DATASET_V1 = zeros(COUNTER, 6);

TRAINING_DATASET_V2 = zeros(COUNTER, 6);


% Common simulation settings.
t0 = 122.83;
% time = 0:stepsize:400;
time = t0:stepsize:t0+400;

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);


% Parallel loop: Each iteration runs an independent simulation.

parfor i = 1:COUNTER

    % Randomize parameters for this simulation.

    L = L_l;% + rand() * (L_u - L_l);

    Amp = A_l;% + rand() * (A_u - A_l);



    % Initialize variables.

%     IC = [-40, 0.01];
%     IC = [-27.5441,0.128307];
    IC = [-20, 0.01];
%     t0 = 271.09;
%     t0 = 122.83;
%     t0 = 20;

% % %     item = 1;
% % % 
    input = zeros(1,2);
% % % 
% % %     input(1) = t0;
% % % 
% % %     %--- Stage 1: Drive system to steady state ---
% % % 
% % %     [~, Y] = ode15s(@FN_MorrisLecar, time, IC, options, item, input, ...
% % %         gK, gCa, gL, EK, EL, ECa, C, V1, V2, V3, V4, phi, Amp, L,I0);
% % % 
% % %     IC = [Y(end, 1), Y(end, 2)];



    %--- Stage 2: Apply first pulse and detect peaks ---

    item = 2;

    input(1) = t0;

    [~, Y] = ode15s(@FN_MorrisLecar, time, IC, options, item, input, gK, gCa, gL, EK, EL, ECa, C, V1, V2, V3, V4, phi, Amp, L, I0);



    [peaks, locs] = findpeaks(Y(:,1), 'MinPeakProminence', 0.05);

    [out, idx] = sort(peaks);

%     spike_location = find(peaks>10); % Spike threshold is defined as 10(mv)


    % Initialize peak-related variables.

    t1 = nan; p1 = nan; t2 = nan; p2 = nan; 
%     t3 = nan; p3 = nan; % Voltage responses from the second pulse
%     I1=nan; % Only needed when we are moving the second pulse

    L1 = nan; L2 = nan;
%     tnextmax = nan; t3max = nan; p3max = nan; l3max = nan; 
%     related to the voltage response corresponding to the second pulse

    if length(locs) > 1

        t1 = time(locs(idx(end)));    % Time of largest peak.

        p1 = out(end);                % Largest peak value.

        I1 = idx(end);                % Index corresponding to largest peak.

        L1 = t1 - t0;



        t2 = time(locs(idx(end-1)));  % Time of second largest peak.

        p2 = out(end-1);

        L2 = t2 - t1;

    end

%% Plotting square forcing
% % % % % % %     if item == 2
% % % % % % % 
% % % % % % %         FORCING1 = I0*ones(1, length(time));
% % % % % % %         L_1 = find(time>=t0); L_2 = find(time<=t0+L);
% % % % % % %         [val,~]=intersect(L_1,L_2);
% % % % % % %         FORCING1(val(1):val(end)) = I0+Amp;
% % % % % % % 
% % % % % % %         figure(1) %first square pulse
% % % % % % %         plot(time,Y(:,1), '-k','linewidth', lw);
% % % % % % %         hold on
% % % % % % %         plot(time(locs), Y(locs), 'or')
% % % % % % %         plot(time, FORCING1 - I0 - 80, '-r', 'linewidth', lw);
% % % % % % % 
% % % % % % %         if isempty(spike_location)
% % % % % % %             text(140, 20, 'subthreshold', 'FontSize', 18)
% % % % % % %         else
% % % % % % %             text(140, 20, 'spike', 'FontSize', 18)
% % % % % % %         end
% % % % % % % 
% % % % % % %         xlabel('t'); ylabel('V(t)');
% % % % % % %         title(sprintf('Amp = %.2f -  DC = %.1f', Amp, php), 'interpreter', 'latex');
% % % % % % %         set(gca, 'FontSize', 18)
% % % % % % %         %xlim([time(1) time(end)/5])
% % % % % % %         %ylim([-90 40])
% % % % % % %         hold off
% % % % % % % 
% % % % % % %         pause(0.1)
% % % % % % % 
% % % % % % %     end
%%
    %--- Stage 3: Adjust the second pulse timing ---

% % % % % % %     item = 3;
% % % % % % % 
% % % % % % %     tnext = t0 + L + L/100;
% % % % % % % 
% % % % % % %     condition = ~(isnan(p1) && isnan(p2));
% % % % % % % 
% % % % % % % 
% % % % % % % 
% % % % % % %     maxOuterIter = 1000;
% % % % % % % 
% % % % % % %     outerIter = 0;
% % % % % % % 
% % % % % % %     while condition && outerIter < maxOuterIter
% % % % % % % 
% % % % % % %         outerIter = outerIter + 1;
% % % % % % % 
% % % % % % %         input(1) = t0;
% % % % % % % 
% % % % % % %         input(2) = tnext;
% % % % % % % 
% % % % % % %         [~, Y] = ode15s(@FN_MorrisLecar, time, IC, options, item, input, ...
% % % % % % %             gK, gCa, gL, EK, EL, ECa, C, V1, V2, V3, V4, phi, Amp, L, I0);
% % % % % % % 
% % % % % % %         [peaks, locs] = findpeaks(Y(:,1), 'MinPeakProminence', 0.05);
% % % % % % % 
% % % % % % %         [out, idx] = sort(peaks);
% % % % % % % 
% % % % % % % 
% % % % % % % 
% % % % % % %         if length(locs) > 1
% % % % % % % 
% % % % % % %             if idx(end) == I1
% % % % % % % 
% % % % % % %                 t3 = time(locs(idx(end-1)));
% % % % % % % 
% % % % % % %                 p3 = out(end-1);
% % % % % % % 
% % % % % % %                 L3 = t3 - t1;
% % % % % % % 
% % % % % % %             else
% % % % % % % 
% % % % % % %                 t3 = time(locs(idx(end)));
% % % % % % % 
% % % % % % %                 p3 = out(end);
% % % % % % % 
% % % % % % %                 L3 = t3 - t1;
% % % % % % % 
% % % % % % %             end
% % % % % % % 
% % % % % % %         end
% % % % % % % 
% % % % % % % 
% % % % % % % 
% % % % % % %         if ~isnan(p2) && ~isnan(p3)
% % % % % % % 
% % % % % % %             maxInnerIter = 1000;
% % % % % % % 
% % % % % % %             innerIter = 0;
% % % % % % % 
% % % % % % %             while p1 > p3 && innerIter < maxInnerIter
% % % % % % % 
% % % % % % %                 innerIter = innerIter + 1;
% % % % % % % 
% % % % % % %                 tnext = tnext + stepsize;
% % % % % % % 
% % % % % % %                 input(1) = t0;
% % % % % % % 
% % % % % % %                 input(2) = tnext;
% % % % % % % 
% % % % % % %                 [~, Y] = ode15s(@FN_MorrisLecar, time, IC, options, item, input, ...
% % % % % % %                     gK, gCa, gL, EK, EL, ECa, C, V1, V2, V3, V4, phi, Amp, L,I0);
% % % % % % % 
% % % % % % %                 [peaks, locs] = findpeaks(Y(:,1), 'MinPeakProminence', 0.05);
% % % % % % % 
% % % % % % %                 [out, idx] = sort(peaks);
% % % % % % % 
% % % % % % %                 if length(locs) > 1
% % % % % % % 
% % % % % % %                     if idx(end) == I1
% % % % % % % 
% % % % % % %                         t3 = time(locs(idx(end-1)));
% % % % % % % 
% % % % % % %                         p3 = out(end-1);
% % % % % % % 
% % % % % % %                         L3 = t3 - t1;
% % % % % % % 
% % % % % % %                     else
% % % % % % % 
% % % % % % %                         t3 = time(locs(idx(end)));
% % % % % % % 
% % % % % % %                         p3 = out(end);
% % % % % % % 
% % % % % % %                         L3 = t3 - t1;
% % % % % % % 
% % % % % % %                     end
% % % % % % % 
% % % % % % %                 end
% % % % % % % 
% % % % % % %                 %% plotting
% % % % % % % % % %                 FORCING1 = I0*ones(1, length(time));
% % % % % % % % % %                 L_1 = find(time>=t0); L_2 = find(time<=t0+L);
% % % % % % % % % %                 [val,~]=intersect(L_1,L_2);
% % % % % % % % % %                 FORCING1(val(1):val(end)) = I0+Amp;
% % % % % % % % % % 
% % % % % % % % % % 
% % % % % % % % % %                 FORCING2 = I0*ones(1, length(time));
% % % % % % % % % %                 L_1 = find(time>=tnext); L_2 = find(time<=tnext+L+L/100);
% % % % % % % % % %                 [val,~]=intersect(L_1,L_2);
% % % % % % % % % %                 FORCING2(val(1):val(end)) = I0+Amp; %error:Index exceeds the number of array elements. Index must not exceed 0.
% % % % % % % % % % %                 length(time)
% % % % % % % % % %                 val(end);
% % % % % % % % % % 
% % % % % % % % % %                 FORCING3 = FORCING1 + FORCING2;
% % % % % % % % % % 
% % % % % % % % % %                 figure(2)
% % % % % % % % % %                 plot(time, FORCING3 - 2*I0 - 80, '-r', 'linewidth', lw);
% % % % % % % % % %                 hold on
% % % % % % % % % % 
% % % % % % % % % % 
% % % % % % % % % %                 plot(time,Y(:,1), '-k','linewidth', lw);
% % % % % % % % % %                 hold on
% % % % % % % % % %                 plot(time(locs), Y(locs), 'or')
% % % % % % % % % % 
% % % % % % % % % %                 if isempty(spike_location)
% % % % % % % % % %                     text(140, 20, 'subthreshold', 'FontSize', 18)
% % % % % % % % % %                 else
% % % % % % % % % %                     text(140, 20, 'spike', 'FontSize', 18)
% % % % % % % % % %                 end
% % % % % % % % % % 
% % % % % % % % % %                 xlabel('t'); ylabel('V(t)');
% % % % % % % % % %                 title(sprintf('P1 = %.3f -  P2 = %.3f - P3 = %.3f', p1, p2, p3), 'interpreter', 'latex');
% % % % % % % % % %                 set(gca, 'FontSize', 18)
% % % % % % % % % %                 %xlim([time(1) time(end)/5])
% % % % % % % % % %                 %ylim([-90 -40])
% % % % % % % % % %                 hold off
% % % % % % % % % % 
% % % % % % % % % %                 pause(0.1)
% % % % % % % 
% % % % % % % 
% % % % % % %             end
% % % % % % % 
% % % % % % %             condition = false; % Exit outer loop if inner loop succeeds.
% % % % % % % %             if innerIter >= maxInnerIter
% % % % % % % %                 warning('Inner loop reached maximum iterations; breaking out.');
% % % % % % % %             end
% % % % % % % %             condition = 0;
% % % % % % % 
% % % % % % %         else
% % % % % % % 
% % % % % % %             break;
% % % % % % % 
% % % % % % %         end
% % % % % % % 
% % % % % % %     end
% % % % % % % 
% % % % % % % 
% % % % % % % 
% % % % % % %     if outerIter >= maxOuterIter
% % % % % % % 
% % % % % % %         warning('Outer loop reached maximum iterations; proceeding with current values.');
% % % % % % % 
% % % % % % %     end
% % % % % % % 
% % % % % % % 
% % % % % % % 
% % % % % % %     %--- Stage 4: Fine-tune the second pulse timing ---
% % % % % % % 
% % % % % % %     pmax = p3;
% % % % % % % 
% % % % % % %     counter_inside = 0;
% % % % % % % 
% % % % % % %     for kx = 1:1000
% % % % % % % 
% % % % % % %         tnext = tnext + stepsize;
% % % % % % % 
% % % % % % %         input(1) = t0;
% % % % % % % 
% % % % % % %         input(2) = tnext;
% % % % % % % 
% % % % % % %         [~, Y] = ode15s(@FN_MorrisLecar, time, IC, options, item, input, ...
% % % % % % %             gK, gCa, gL, EK, EL, ECa, C, V1, V2, V3, V4, phi, Amp, L,I0);
% % % % % % % 
% % % % % % %         [peaks, locs] = findpeaks(Y(:,1), 'MinPeakProminence', 0.05);
% % % % % % % 
% % % % % % %         [out, idx] = sort(peaks);
% % % % % % % 
% % % % % % %         if length(locs) > 1
% % % % % % % 
% % % % % % %             if idx(end) == I1
% % % % % % % 
% % % % % % %                 t3 = time(locs(idx(end-1)));
% % % % % % % 
% % % % % % %                 p3 = out(end-1);
% % % % % % % 
% % % % % % %                 L3 = t3 - t1;
% % % % % % % 
% % % % % % %             else
% % % % % % % 
% % % % % % %                 t3 = time(locs(idx(end)));
% % % % % % % 
% % % % % % %                 p3 = out(end);
% % % % % % % 
% % % % % % %                 L3 = t3 - t1;
% % % % % % % 
% % % % % % %             end
% % % % % % % 
% % % % % % %         end
% % % % % % % 
% % % % % % % 
% % % % % % % 
% % % % % % %         if p3 > pmax
% % % % % % % 
% % % % % % %             t3max = t3;
% % % % % % % 
% % % % % % %             p3max = p3;
% % % % % % % 
% % % % % % %             pmax = p3;
% % % % % % % 
% % % % % % %             l3max = L3;
% % % % % % % 
% % % % % % %             tnextmax = tnext;
% % % % % % % 
% % % % % % %             counter_inside = 0;
% % % % % % % 
% % % % % % %         else
% % % % % % % 
% % % % % % %             counter_inside = counter_inside + 1;
% % % % % % % 
% % % % % % %         end
% % % % % % % 
% % % % % % % %% Plotting
% % % % % % % % % %         FORCING1 = I0*ones(1, length(time));
% % % % % % % % % %         L_1 = find(time>=t0); L_2 = find(time<=t0+L);
% % % % % % % % % %         [val,~]=intersect(L_1,L_2);
% % % % % % % % % %         FORCING1(val(1):val(end)) = I0+Amp;
% % % % % % % % % % 
% % % % % % % % % % 
% % % % % % % % % %         FORCING2 = I0*ones(1, length(time));
% % % % % % % % % %         L_1 = find(time>=tnext); L_2 = find(time<=tnext+L+L/100);
% % % % % % % % % %         [val,~]=intersect(L_1,L_2);
% % % % % % % % % %         FORCING2(val(1):val(end)) = I0+Amp;
% % % % % % % % % % 
% % % % % % % % % %         FORCING3 = FORCING1 + FORCING2;
% % % % % % % % % % 
% % % % % % % % % %         figure(3)
% % % % % % % % % %         plot(time, FORCING3 - 2*I0 - 80, '-r', 'linewidth', lw);
% % % % % % % % % %         hold on
% % % % % % % % % % 
% % % % % % % % % % 
% % % % % % % % % %         plot(time,Y(:,1), '-k','linewidth', lw);
% % % % % % % % % %         hold on
% % % % % % % % % %         plot(time(locs), Y(locs), 'or')
% % % % % % % % % % 
% % % % % % % % % %         if isempty(spike_location)
% % % % % % % % % %             text(140, 20, 'subthreshold', 'FontSize', 18)
% % % % % % % % % %         else
% % % % % % % % % %             text(140, 20, 'spike', 'FontSize', 18)
% % % % % % % % % %         end
% % % % % % % % % % 
% % % % % % % % % %         xlabel('t'); ylabel('V(t)');
% % % % % % % % % %         title(sprintf('P1 = %.3f -  P2 = %.3f - P3 = %.3f', p1, p2, p3), 'interpreter', 'latex');
% % % % % % % % % %         set(gca, 'FontSize', 18)
% % % % % % % % % %         %xlim([time(1) time(end)/5])
% % % % % % % % % %         %ylim([-90 -40])
% % % % % % % % % %         hold off
% % % % % % % % % % 
% % % % % % % % % % 
% % % % % % % % % %         pause(0.1)
% % % % % % %         %%
% % % % % % % 
% % % % % % % 
% % % % % % %         if counter_inside > 20
% % % % % % % 
% % % % % % %             break;
% % % % % % % 
% % % % % % %         end
% % % % % % % 
% % % % % % %     end



    % Store results from this iteration.

% % % % % %     TRAINING_DATASET_V1(i, :) = [Amp, L, t1, p1, t2, p2, tnextmax, t3max, p3max];
% % % % % % 
% % % % % %     TRAINING_DATASET_V2(i, :) = [Amp, L, L1, p1, L2, p2, tnextmax, l3max, p3max];

    TRAINING_DATASET_V1(i, :) = [Amp, L, t1, p1, t2, p2];

    TRAINING_DATASET_V2(i, :) = [Amp, L, L1, p1, L2, p2];

end


% Save the aggregated results to .mat files after the parfor loop.

save(sprintf('TRAINING_DATASET_V1_%d.mat', seed), 'TRAINING_DATASET_V1');

save(sprintf('TRAINING_DATASET_V2_%d.mat', seed), 'TRAINING_DATASET_V2');

% disp('TRAINING_DATASET_V1 =');
% disp(TRAINING_DATASET_V1);

disp('TRAINING_DATASET_V2 =');
disp(TRAINING_DATASET_V2);
end