% This code is supplementary material of the publication
%
% Denk, F., Kollmeier, B and Ewert, S.D.
% "Removing reflections in semianechoic impulse responses by 
%  frequency-dependent truncation "
% Journal of the Audio Engineering Society, vol 66, No 4, 2018
% DOI: https://doi.org/10.17743/jaes.2018.0002
%
% Please cite the paper when using the code in this repository.
%
% This software comes free with the hope that it is useful, but without any
% warranty. 
% It is published under the terms of the GNU General Public
% License v3.0. You are free to use, modify and redistribute the code,
% provided the original source is attributed and further distribution is
% made under the same license.
%
%% Copyright (c) Florian Denk, 2018
% Email: florian.denk@uni-oldenburg.de
% Department of Medical Physics and Acoustics, University of Oldenburg

%%
% This scripts shows the processing steps of Experiment II and replicates
% Figures 4 and 5
%
% Florian Denk, 21.12.2017

clear all
close all
clc

%% General Parameters


n_fft = 8128; % FFT length for Plotting

% Colors and line marks
black =        [0 0 0];
grey1 = .45 .* [1 1 1];
grey2 = .75 .* [1 1 1];

cl_an   = grey2;
cl_semi = black;
cl_bt   = grey2;
cl_stft = grey1;
cl_dft  = grey1;
cl_sm   = black;

mr_an   = '-';
mr_semi = '-';
mr_bt   = '-.';
mr_stft = '-';
mr_dft  = '-.';
mr_sm   = '-.';

lw_sm = 1;
lw_b = 2;


%% Load Measured IRs

S_Semi = load('Data\Semi_Anechoic.mat');
S_Anec = load('Data\Anechoic.mat');

v_ir_semi = S_Semi.v_ir;
v_ir_anec = S_Anec.v_ir;

% Level Adjustments for plots
v_ir_semi = v_ir_semi .* db2mag(15);
v_ir_anec = v_ir_anec .* db2mag(15);

% Read sampling rate
srate = S_Semi.srate;

% Frequency Vector
v_frq = linspace(0, srate/2, n_fft/2 + 1)';


%% Processing Parameters

% STFT parameters for FDT
n_block =  180; % 3.33 ms - 300 Hz frequncy resolution
n_shift =   90; % 0.33 ms

% Generate STFT frequency vector
v_frq_stft = linspace(0, srate/2, n_block/2+1)';

% Shortest truncation length in FDT = BT truncation length
T_min = 0.003;    % in seconds
T_tr_BT = 0.003;  

% Frequency dependent truncation length
v_trwin_lengths = zeros(size(v_frq_stft));
v_trwin_lengths(1)     = 0.200; %   0 - 150 Hz
v_trwin_lengths(2)     = 0.016; % 150 - 450 Hz
v_trwin_lengths(3)     = 0.006; % 450 - 600 Hz
v_trwin_lengths(4)     = 0.004; % 600 - 750 Hz
v_trwin_lengths(5:end) = T_min; % > 750 Hz

% Read delay of impulse response
[~,n_delay_ir] = max(abs(hilbert(v_ir_semi)));
if mod(n_delay_ir,2)
    n_delay_ir = n_delay_ir +1;
end

%% Processing

% Broadband truncation
v_ir_bt = BBTruncate(v_ir_semi, 1, round( srate*T_tr_BT) + n_delay_ir , 0, 24);

% FDT-STFT (This is what we recommend)
[v_ir_fdt, S_FDTProcSteps] = ...
                    FDTruncate(v_ir_semi, ...
                               v_trwin_lengths,...
                               srate, ...
                               n_block, ...
                               n_shift);

% bring all to same length
if length(v_ir_anec) <= n_fft
   v_ir_anec = [v_ir_anec;   zeros(n_fft - size(v_ir_anec ,1), 1)];
else
   v_ir_anec = BBTruncate(v_ir_anec,1,n_fft,0,48);
end
if length(v_ir_semi) <= n_fft
   v_ir_semi = [v_ir_semi;   zeros(n_fft - size(v_ir_semi ,1), 1)];
else
   v_ir_semi = BBTruncate(v_ir_semi,1,n_fft,0,48);
end
if length(v_ir_bt) <= n_fft
   v_ir_bt = [v_ir_bt;   zeros(n_fft - size(v_ir_bt ,1), 1)];
else
   v_ir_bt = BBTruncate(v_ir_bt,1,n_fft,0,48);
end
if length(v_ir_fdt) <= n_fft
   v_ir_fdt = [v_ir_fdt;   zeros(n_fft - size(v_ir_fdt ,1), 1)];
else
   v_ir_fdt = BBTruncate(v_ir_fdt,1,n_fft,0,48);
end


% Time vector
v_t_ir = (0:1/srate: 1/srate * (n_fft-1)) - ((n_delay_ir-1)/srate);

%% Figure 4


% Write IRs together and set plot parameters
M_IR(:,1) = v_ir_semi;
M_IR(:,2) = v_ir_anec;

v_flim = [20 20000];
v_tlim = [-0.002 0.021];
v_fticks = [30 62 125 250 500 1000 2000 4000 8000 16000];
v_speclim = [-23 5];

plot_wins = 1;
plot_IRlog = 1;
M_colors = [cl_semi;cl_an];

% Summary Plot
plot_FDTsummary(v_t_ir,...
                M_IR,...
                S_FDTProcSteps,...
                n_fft,...
                {'Semi-Anech.', 'Anech.'},...
                v_flim,...
                v_tlim,...
                v_fticks,...
                v_speclim,...
                plot_wins,...
                plot_IRlog,...
                M_colors);
drawnow;
set(gcf, 'Name', 'Fig. 4')



%% Figure 5

% IR Comparison
figure('position',[50 50 800 800], 'Name', 'Fig. 5')

subplot(211)
plot(v_t_ir*1000, mag2db(abs(v_ir_semi)), mr_semi, 'Color', cl_semi, 'LineWidth', lw_b); hold on
plot(v_t_ir*1000, mag2db(abs(v_ir_anec)), mr_an,   'Color', cl_an  , 'LineWidth', lw_sm);   hold on
plot(v_t_ir*1000, mag2db(abs(v_ir_fdt )), mr_stft, 'Color', cl_stft, 'LineWidth', lw_sm); hold on
plot(v_t_ir*1000, mag2db(abs(v_ir_bt  )), mr_bt,   'Color', cl_bt, 'LineWidth', lw_sm);   hold on
grid on
box on
xlim([-1 20])
ylim([-100 0])
xlabel('t [ms]')
ylabel('log(|h|) [dB]')
set(gca, 'Ytick',[-90:10:-10],...
         'YTickLabel',strread(num2str([-90:10:10]), '%s')'  );
legend('Semi-Anech.', 'Anech.', 'FDT', 'BT', 'Location', 'NorthEast',...
    'Orientation', 'Horizontal');

% Response Comparison
subplot(212)
semilogx(v_frq, mag2db(abs( fftR( v_ir_semi) )), mr_semi, 'Color', cl_semi, 'LineWidth', lw_b); hold on
semilogx(v_frq, mag2db(abs( fftR( v_ir_anec) )), mr_an,   'Color', cl_an  , 'LineWidth', lw_sm); 
semilogx(v_frq, mag2db(abs( fftR( v_ir_fdt)  )), mr_stft, 'Color', cl_stft, 'LineWidth', lw_sm); 
semilogx(v_frq, mag2db(abs( fftR( v_ir_bt)   )), mr_bt,   'Color', cl_bt  , 'LineWidth', lw_sm); 
grid on
box on
xlim([30 24000])
ylim([-25 5])
xlabel('f [kHz]')
ylabel('|H| [dB]')
set(gca, 'Ytick',[-25:5:5],...
         'YtickLabel', {'' '-20' '' '-10' '' '0' ''},...
         'XTick', [0.1 0.4 1 2 4 6 8 12 16] * 1000,...
         'XTickLabel', {'0.1', '0.4', '1', '2', '4', '6', '8', '12', '16'});
    

legend('Semi-Anech.', 'Anech.', 'FDT', 'BT', 'Location', 'NorthEast',...
    'Orientation', 'Horizontal');

