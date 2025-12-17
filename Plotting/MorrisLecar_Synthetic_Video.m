function [TRAINING_DATASET_V1, TRAINING_DATASET_V2] = MorrisLecar_Synthetic_Video(seed, COUNTER)

% ML_Data_Driven5_parallel - Generate data based on the Morris–Lecar model

% using MATLAB's parallel computing capabilities.

%

% Usage:

%   [DS1, DS2] = ML_Data_Driven5_parallel(seed, COUNTER)

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

A_l = 15;
L_l = 40;
I0=80;

% Use provided COUNTER or default to 1 if not given

if nargin < 2

    COUNTER = 1;

end



% Set the random seed for reproducibility.

% rng(seed, 'twister');

% Preallocate output arrays.

TRAINING_DATASET_V1 = zeros(COUNTER, 9);

TRAINING_DATASET_V2 = zeros(COUNTER, 9);



% Common simulation settings.

time = 0:stepsize:400;

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);


% Parallel loop: Each iteration runs an independent simulation.

% parfor i = 1:COUNTER
for i = 1:COUNTER



    % Randomize parameters for this simulation.

    L = L_l;% + rand() * (L_u - L_l);

    Amp = A_l;% + rand() * (A_u - A_l);



    % Initialize variables.

    IC = [-20, 0.01];
%     IC = [-27, 0.1];

    item = 1;

    t0 = 20;
%     t0 = 131.4;
%     t0 = 92;


    input = zeros(1,2);

    input(1) = t0;



    %--- Stage 1: Drive system to steady state ---

    [~, Y] = ode15s(@FN_MorrisLecar, time, IC, options, item, input, ...
        gK, gCa, gL, EK, EL, ECa, C, V1, V2, V3, V4, phi, Amp, L,I0);

    IC = [Y(end, 1), Y(end, 2)];



    %--- Stage 2: Apply first pulse and detect peaks ---

    item = 2;

    input(1) = t0;

    [~, Y] = ode15s(@FN_MorrisLecar, time, IC, options, item, input, gK, gCa, gL, EK, EL, ECa, C, V1, V2, V3, V4, phi, Amp, L, I0);



    [peaks, locs] = findpeaks(Y(:,1), 'MinPeakProminence', 0.05);

    [out, idx] = sort(peaks);

    spike_location = find(peaks>10); % Spike threshold is defined as 10(mv)


    % Initialize peak-related variables.

    t1 = nan; p1 = nan; t2 = nan; p2 = nan; t3 = nan; p3 = nan;
    I1=nan;

    tnextmax = nan; t3max = nan; p3max = nan; L1 = nan; L2 = nan; l3max = nan;



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
    if item == 2

        FORCING1 = I0*ones(1, length(time));
        L_1 = find(time>=t0); L_2 = find(time<=t0+L);
        [val,~]=intersect(L_1,L_2);
        FORCING1(val(1):val(end)) = I0+Amp;

% % % %         figure(1) %first square pulse
% % % %         plot(time,Y(:,1), '-k','linewidth', lw);
% % % %         hold on
% % % %         plot(time(locs), Y(locs), 'or')
% % % %         plot(time, FORCING1 - I0 - 80, '-r', 'linewidth', lw);
% % % % 
% % % %         if isempty(spike_location)
% % % %             text(140, 20, 'subthreshold', 'FontSize', 18)
% % % %         else
% % % %             text(140, 20, 'spike', 'FontSize', 18)
% % % %         end
% % % % 
% % % %         xlabel('t'); ylabel('V(t)');
% % % %         title(sprintf('Amp = %.2f -  DC = %.1f', Amp, php), 'interpreter', 'latex');
% % % %         set(gca, 'FontSize', 18)
% % % %         %xlim([time(1) time(end)/5])
% % % %         %ylim([-90 40])
% % % %         hold off
% % % % 
% % % %         pause(0.1)

    end
%%
    %--- Stage 3: Adjust the second pulse timing ---

    item = 3;

    tnext = t0 + L + L/100;

    condition = ~(isnan(p1) && isnan(p2));



    maxOuterIter = 1000;

    outerIter = 0;

    %% === VIDEO SETUP ===
    v = VideoWriter('spike_animation.mp4','MPEG-4'); % Save as MP4
    v.FrameRate = 10; % adjust speed (frames per second)
    open(v);

    while condition && outerIter < maxOuterIter

        outerIter = outerIter + 1;

        input(1) = t0;

        input(2) = tnext;

        [~, Y] = ode15s(@FN_MorrisLecar, time, IC, options, item, input, ...
            gK, gCa, gL, EK, EL, ECa, C, V1, V2, V3, V4, phi, Amp, L, I0);

        [peaks, locs] = findpeaks(Y(:,1), 'MinPeakProminence', 0.05);

        [out, idx] = sort(peaks);



        if length(locs) > 1

            if idx(end) == I1

                t3 = time(locs(idx(end-1)));

                p3 = out(end-1);

                L3 = t3 - t1;

            else

                t3 = time(locs(idx(end)));

                p3 = out(end);

                L3 = t3 - t1;

            end

        end



        if ~isnan(p2) && ~isnan(p3)

            maxInnerIter = 1000;

            innerIter = 0;

            while p1 > p3 && innerIter < maxInnerIter

                innerIter = innerIter + 1;

                tnext = tnext + stepsize;

                input(1) = t0;

                input(2) = tnext;

                [~, Y] = ode15s(@FN_MorrisLecar, time, IC, options, item, input, ...
                    gK, gCa, gL, EK, EL, ECa, C, V1, V2, V3, V4, phi, Amp, L,I0);

                [peaks, locs] = findpeaks(Y(:,1), 'MinPeakProminence', 0.05);

                [out, idx] = sort(peaks);

                if length(locs) > 1

                    if idx(end) == I1

                        t3 = time(locs(idx(end-1)));

                        p3 = out(end-1);

                        L3 = t3 - t1;

                    else

                        t3 = time(locs(idx(end)));

                        p3 = out(end);

                        L3 = t3 - t1;

                    end

                end

                %% plotting
                FORCING1 = I0*ones(1, length(time));
                L_1 = find(time>=t0); L_2 = find(time<=t0+L);
                [val,~]=intersect(L_1,L_2);
                FORCING1(val(1):val(end)) = I0+Amp;


                FORCING2 = I0*ones(1, length(time));
                L_1 = find(time>=tnext); L_2 = find(time<=tnext+L+L/100);
                [val,~]=intersect(L_1,L_2);
                FORCING2(val(1):val(end)) = I0+Amp; %error:Index exceeds the number of array elements. Index must not exceed 0.
%                 length(time)
                val(end);

                FORCING3 = FORCING1 + FORCING2;

% % %                 figure(2)
% % %                 plot(time, FORCING3 - 2*I0 - 80, '-r', 'linewidth', lw);
% % %                 hold on
% % % 
% % % 
% % %                 plot(time,Y(:,1), '-k','linewidth', lw);
% % %                 hold on
% % %                 plot(time(locs), Y(locs), 'or')
% % % 
% % %                 if isempty(spike_location)
% % %                     text(140, 20, 'subthreshold', 'FontSize', 18)
% % %                 else
% % %                     text(140, 20, 'spike', 'FontSize', 18)
% % %                 end
% % % 
% % %                 xlabel('t'); ylabel('V(t)');
% % %                 title(sprintf('P1 = %.3f -  P2 = %.3f - P3 = %.3f', p1, p2, p3), 'interpreter', 'latex');
% % %                 set(gca, 'FontSize', 18)
% % %                 %xlim([time(1) time(end)/5])
% % %                 %ylim([-90 -40])
% % %                 hold off
% % % 
% % %                 pause(0.1)

                
                % For video (Same as figure(2). Hence figure (2) is
                % commented
                figure(4)
                subplot(4,1,[1 3])
                plot(time,Y(:,1), '-b','linewidth', lw);
                hold on
                plot(time(locs), Y(locs), 'or')
                ylabel('$V$',Interpreter='latex')
                title(sprintf('P1 = %.3f -  P2 = %.3f - P3 = %.3f', p1, p2, p3), 'interpreter', 'latex');
                set(groot,'defaultAxesTickLabelInterpreter','latex');
                set(groot,'defaulttextinterpreter','latex');
                set(groot,'defaultLegendInterpreter','latex');
                set(gca,'TickLabelInterpreter','latex','fontsize',18);
                xticklabels([]); % Remove x-axis labels in this subplot
% set(gca, 'FontSize', 18)
                %xlim([time(1) time(end)/5])
                ylim([-60 50])
                hold off

                subplot(4,1,4)
                plot(time, FORCING3-I0, '-r', 'linewidth', lw);
%                 hold on
                set(groot,'defaultAxesTickLabelInterpreter','latex');
                set(groot,'defaulttextinterpreter','latex');
                set(groot,'defaultLegendInterpreter','latex');
                set(gca,'TickLabelInterpreter','latex','fontsize',18);
                ylabel('$I_{app}$',Interpreter='latex')
                xlabel('$t$',Interpreter='latex')
                yticks(80:15:95)  % ticks every 25 mV

                % === CAPTURE FRAME FOR VIDEO ===
                frame = getframe(gcf);
                writeVideo(v, frame);
            end
           
            % === CLOSE VIDEO FILE ===
            close(v);
            disp('Video saved as spike_animation.mp4');

            condition = false; % Exit outer loop if inner loop succeeds.
%             if innerIter >= maxInnerIter
%                 warning('Inner loop reached maximum iterations; breaking out.');
%             end
%             condition = 0;

        else

            break;

        end

    end



    if outerIter >= maxOuterIter

        warning('Outer loop reached maximum iterations; proceeding with current values.');

    end



    %--- Stage 4: Fine-tune the second pulse timing ---

    pmax = p3;

    counter_inside = 0;

    for kx = 1:1000

        tnext = tnext + stepsize;

        input(1) = t0;

        input(2) = tnext;

        [~, Y] = ode15s(@FN_MorrisLecar, time, IC, options, item, input, ...
            gK, gCa, gL, EK, EL, ECa, C, V1, V2, V3, V4, phi, Amp, L,I0);

        [peaks, locs] = findpeaks(Y(:,1), 'MinPeakProminence', 0.05);

        [out, idx] = sort(peaks);

        if length(locs) > 1

            if idx(end) == I1

                t3 = time(locs(idx(end-1)));

                p3 = out(end-1);

                L3 = t3 - t1;

            else

                t3 = time(locs(idx(end)));

                p3 = out(end);

                L3 = t3 - t1;

            end

        end



        if p3 > pmax

            t3max = t3;

            p3max = p3;

            pmax = p3;

            l3max = L3;

            tnextmax = tnext;

            counter_inside = 0;

        else

            counter_inside = counter_inside + 1;

        end

%% Plotting
        FORCING1 = I0*ones(1, length(time));
        L_1 = find(time>=t0); L_2 = find(time<=t0+L);
        [val,~]=intersect(L_1,L_2);
        FORCING1(val(1):val(end)) = I0+Amp;


        FORCING2 = I0*ones(1, length(time));
        L_1 = find(time>=tnext); L_2 = find(time<=tnext+L+L/100);
        [val,~]=intersect(L_1,L_2);
        FORCING2(val(1):val(end)) = I0+Amp;

        FORCING3 = FORCING1 + FORCING2;
% % % % 
% % % %         figure(3)
% % % %         plot(time, FORCING3 - 2*I0 - 80, '-r', 'linewidth', lw);
% % % %         hold on
% % % % 
% % % % 
% % % %         plot(time,Y(:,1), '-k','linewidth', lw);
% % % %         hold on
% % % %         plot(time(locs), Y(locs), 'or')
% % % % 
% % % %         if isempty(spike_location)
% % % %             text(140, 20, 'subthreshold', 'FontSize', 18)
% % % %         else
% % % %             text(140, 20, 'spike', 'FontSize', 18)
% % % %         end
% % % % 
% % % %         xlabel('t'); ylabel('V(t)');
% % % %         title(sprintf('P1 = %.3f -  P2 = %.3f - P3 = %.3f', p1, p2, p3), 'interpreter', 'latex');
% % % %         set(gca, 'FontSize', 18)
% % % %         %xlim([time(1) time(end)/5])
% % % %         %ylim([-90 -40])
% % % %         hold off
% % % % 
% % % % 
% % % %         pause(0.1)
        %%


        if counter_inside > 20

            break;

        end

    end



    % Store results from this iteration.

    TRAINING_DATASET_V1(i, :) = [Amp, L, t1, p1, t2, p2, tnextmax, t3max, p3max];

    TRAINING_DATASET_V2(i, :) = [Amp, L, L1, p1, L2, p2, tnextmax, l3max, p3max];

end


% Save the aggregated results to .mat files after the parfor loop.

save(sprintf('TRAINING_DATASET_V1_%d.mat', seed), 'TRAINING_DATASET_V1');

save(sprintf('TRAINING_DATASET_V2_%d.mat', seed), 'TRAINING_DATASET_V2');

% disp('TRAINING_DATASET_V1 =');
% disp(TRAINING_DATASET_V1);

disp('TRAINING_DATASET_V2 =');
disp(TRAINING_DATASET_V2);
end
