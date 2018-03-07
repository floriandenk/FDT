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
% This scripts shows the processing steps of Experiment I and replicates
% Figures 1 and 2
%
% Florian Denk, 21.12.2017

clear all
close all
clc

addpath(genpath(pwd))


%% General Parameters

srate = 48000; % Sampling rate

n_fft = 4096; % FFT length for Plotting

% Colors and lines
black =       [0 0 0];
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

%% Construct Dummy Loudspeaker Response

% Frequency limits of bandpass
f_bp = [65 18000];

% Bandpass IR
[b, a] = butter(6, f_bp/(srate/2), 'bandpass');
v_ir_anec = impz(b, a, 5000); clear a b

% Add initial delay
n_delay_first = 480;
v_ir_anec = [zeros(n_delay_first,1); v_ir_anec];

%% Construct Reflection Filter
f_cut = 1000; % lowest frequency with considerable energy
n_refl = 1024; % length of reflection IR

[b_r, a_r] = butter(4, f_cut/(srate/2), 'high');

% IR of reflextion filter
v_refl = filter(b_r, a_r, v_ir_anec(1:n_refl)); clear b_r a_r


%% Construct semi-anechoic IR

% Delay and gain of reflection
gain_refl    = db2mag(-10);
n_start_refl = round(0.006*srate); % 6 ms delay

% Add reflection to anechoic IR
v_ir_semian = v_ir_anec;
v_ir_semian    (n_start_refl:n_start_refl+n_refl-1) =  ...
    v_ir_semian(n_start_refl:n_start_refl+n_refl-1) + (gain_refl * v_refl);


%% Processing Parameters

% STFT parameters for FDT
n_block =  180; % 3.33 ms - 300 Hz frequncy resolution
n_shift =   90; % 0.33 ms

% Generate STFT frequency vector
v_frq_stft = linspace(0, srate/2, n_block/2+1)';

% Shortest truncation length in FDT = BT truncation length
T_min = 0.0045; % in seconds
T_tr_TD = 0.0045;

% Frequency-dependent truncation times
v_trwin_lengths = zeros(size(v_frq_stft));
v_trwin_lengths(1)     = 0.100; %   0 - 150 Hz
v_trwin_lengths(2)     = 0.030; % 150 - 450 Hz
v_trwin_lengths(3)     = 0.010; % 450 - 600 Hz
v_trwin_lengths(4)     = 0.006; % 600 - 750 Hz
v_trwin_lengths(5:end) = T_min; % > 750 Hz



%% Processing

% Broadband Truncation
v_ir_bt = BBTruncate(v_ir_semian, 1, round( srate*(T_tr_TD + n_delay_first/srate) ), 8, 24);

% FDT-STFT (This is what we recommend)
[v_ir_stft, S_ProcSteps] = ...
                    FDTruncate(v_ir_semian, ...
                               v_trwin_lengths,...
                               srate, ...
                               n_block, ...
                               n_shift);
    
% FDT-DFT
v_ir_dft = FDTruncate_DFT(v_ir_semian, v_trwin_lengths, srate);      


% Complex Spectral Smoothing
% For details on complex smmothing, see Hatziantoniou & Mourjopoulos, 
% JAES  2000

% Parameters
S_cfg.v_frq = linspace(0, 24000, n_fft/2+1)';
S_cfg.s_win = 'rect';
S_cfg.s_fscale = 'FDT_eq'; % Equivalent lendths to FDT

% Expand truncation windows to spectral sampling in spectral smoothing -
% but with same spectral resolution, by nearest neighbout interpolation
S_cfg.v_trwin =  interp1(linspace(0, 24000, length(v_frq_stft))',...
                        v_trwin_lengths,...
                        S_cfg.v_frq, ...
                        'nearest');
% Get smoothing matrix: Smoothing window for each frequency bin
M_sm = csmooth_init_FDTeq(S_cfg); 

% Do smoothing
    % Remove initial delay and convert to frequency domain
    v_ir_sm = fftR(circshift(v_ir_semian, -n_delay_first), n_fft);
    % Apply smoothing matrix and get pack to time domain
    v_ir_sm = ifftR( M_sm * v_ir_sm );
    % Reapply delay
    v_ir_sm = circshift(v_ir_sm, n_delay_first);


% Put all to length of fft
if length(v_ir_anec) <= n_fft
   v_ir_anec = [v_ir_anec;   zeros(n_fft - size(v_ir_anec ,1), 1)];
else
   v_ir_anec = v_ir_anec(1:n_fft);
end
if length(v_ir_semian) <= n_fft
   v_ir_semian = [v_ir_semian;   zeros(n_fft - size(v_ir_semian ,1), 1)];
else
   v_ir_semian = v_ir_semian(1:n_fft);
end
if length(v_ir_bt) <= n_fft
   v_ir_bt = [v_ir_bt;   zeros(n_fft - size(v_ir_bt ,1), 1)];
else
   v_ir_bt = v_ir_bt(1:n_fft);
end
if length(v_ir_dft) <= n_fft
   v_ir_dft = [v_ir_dft;   zeros(n_fft - size(v_ir_dft ,1), 1)];
else
   v_ir_dft = v_ir_dft(1:n_fft);
end
if length(v_ir_stft) <= n_fft
   v_ir_stft = [v_ir_stft;   zeros(n_fft - size(v_ir_stft ,1), 1)];
else
   v_ir_stft = v_ir_stft(1:n_fft);
end
if length(v_ir_sm) <= n_fft
   v_ir_sm = [v_ir_sm;   zeros(n_fft - size(v_ir_sm ,1), 1)];
else
   v_ir_sm = v_ir_sm(1:n_fft);
end


% Time vector
v_t_ir = (0:1/48000: 1/48000 * (n_fft-1)) - (n_delay_first-1)/srate;



%% Plot: Figure 1

% Write IRs together
M_IR(:,1) = v_ir_semian;
M_IR(:,2) = v_ir_anec;
M_IR(:,3) = [zeros(n_start_refl-1,1); v_refl .* gain_refl; zeros(n_fft - n_start_refl+1 - length(v_refl), 1)];

if mod(size(M_IR,1),2)
    M_IR = [M_IR; zeros(1,3)];
end



% Better plot representation
% Now in separate function
v_flim = [20 20000];
v_tlim = [-0.002 0.021];
v_fticks = [30 62 125 250 500 1000 2000 4000 8000 16000];
v_speclim = [-23 5];

plot_wins = 1;
plot_IRlog = 1;
M_colors = [cl_semi;cl_an;cl_bt];

plot_FDTsummary(v_t_ir,...
                M_IR,...
                S_ProcSteps,...
                n_fft,...
                {'Semi-Anech.', 'Anech.','Refl.'},...
                v_flim,...
                v_tlim,...
                v_fticks,...
                v_speclim,...
                plot_wins,...
                plot_IRlog,...
                M_colors);
       
set(gcf, 'name', 'Fig. 1');

%% Plot: Figure 2, left panels

xlim_ir = [-1 20]; % In ms!!

figure('name', 'Fig. 2, left', 'position', [50 50 800 800])

subplot(211)
h_an   = plot(v_t_ir*1000, mag2db(abs(v_ir_anec)),   mr_an   , 'Color', cl_an,    'LineWidth', lw_b ); hold on
h_semi = plot(v_t_ir*1000, mag2db(abs(v_ir_semian)), mr_semi , 'Color', cl_semi,  'LineWidth', lw_sm, 'Visible','Off');
h_sm   = plot(v_t_ir*1000, mag2db(abs(v_ir_sm)),     mr_sm   , 'Color', cl_sm,    'LineWidth', lw_sm);
h_bt   = plot(v_t_ir*1000, mag2db(abs(v_ir_bt)),     mr_bt   , 'Color', cl_bt,    'LineWidth', lw_sm);
h_dft  = plot(v_t_ir*1000, mag2db(abs(v_ir_dft)),    mr_dft  , 'Color', cl_dft  , 'LineWidth', lw_sm);
h_stft = plot(v_t_ir*1000, mag2db(abs(v_ir_stft)),   mr_stft , 'Color', cl_stft , 'LineWidth', lw_sm); 
% h_subb = plot(v_t_ir*1000, v_ir_subb,   '-y'   , 'LineWidth', 2); 

box on
xlim(xlim_ir)
set(gca,'YTick', [-90:10:-10],...
        'YLim', [-100 0]);
grid on

ylabel('log(|h|) [dB]')
title('Processed IRs')

legend([h_an h_semi h_sm h_bt h_dft h_stft], ...
       {'Anech.', 'Semi-Anech.', 'Smooth.', 'BT', 'FDT-DFT', 'FDT-STFT'},...
       'Location', 'EastOutside')

subplot(212)
v_snr_semian = min(200, 10*log10( abs(v_ir_anec.^2) ./ abs(v_ir_anec - v_ir_semian).^2 ));
v_snr_bt     = min(200, 10*log10( abs(v_ir_anec.^2) ./ abs(v_ir_anec - v_ir_bt).^2 ));
v_snr_dft    = min(200, 10*log10( abs(v_ir_anec.^2) ./ abs(v_ir_anec - v_ir_dft).^2 ));
v_snr_stft   = min(200, 10*log10( abs(v_ir_anec.^2) ./ abs(v_ir_anec - v_ir_stft).^2 ));
v_snr_sm     = min(200, 10*log10( abs(v_ir_anec.^2) ./ abs(v_ir_anec - v_ir_sm).^2 ));

% Smooth SNR for plotting
n_sm_snr = 48; % Just for plotting
plot(v_t_ir*1000, smooth( v_snr_semian, n_sm_snr), mr_semi , 'Color', cl_semi,  'LineWidth', lw_sm); hold on
plot(v_t_ir*1000, smooth( v_snr_sm    , n_sm_snr), mr_sm   , 'Color', cl_sm,    'LineWidth', lw_sm);
plot(v_t_ir*1000, smooth( v_snr_bt    , n_sm_snr), mr_bt   , 'Color', cl_bt,    'LineWidth', lw_sm);
plot(v_t_ir*1000, smooth( v_snr_dft   , n_sm_snr), mr_dft  , 'Color', cl_dft  , 'LineWidth', lw_sm); hold on
plot(v_t_ir*1000, smooth( v_snr_stft  , n_sm_snr), mr_stft , 'Color', cl_stft , 'LineWidth', lw_sm); hold on

box on
xlim(xlim_ir)
ylim([-45 110])
set(gca, 'Ytick',[-40:20:100],...
         'YtickLabel',{'-40' '-20' '0' '20' '40' '60' '80' '' })
grid on
xlabel('t [ms]')
ylabel('SNR [dB]')
title('Deviation from Anechoic IR')
drawnow



%% Plot Figure 2, right panels
v_frq_stft = linspace(0, 24000, n_fft/2 + 1)';

xlim_pl = [50 20000];


figure('name', 'Fig 2, right', 'position', [50 50 800 800])
subplot(311)
y_sep = 5;
h_semi = semilogx(v_frq_stft, mag2db(abs(fftR(v_ir_semian)))-0*y_sep, mr_semi  , 'Color', cl_semi,  'LineWidth', lw_sm, 'Visible', 'On'); hold on
h_an   = semilogx(v_frq_stft, mag2db(abs(fftR(v_ir_anec)))  -0*y_sep, mr_an    , 'Color', cl_an,    'LineWidth', lw_b); hold on
h_bt   = semilogx(v_frq_stft, mag2db(abs(fftR(v_ir_bt)))    -1*y_sep, mr_bt    , 'Color', cl_bt,    'LineWidth', lw_sm);
h_dft  = semilogx(v_frq_stft, mag2db(abs(fftR(v_ir_dft)))   -2*y_sep, mr_dft   , 'Color', cl_dft,   'LineWidth', lw_sm); 
h_stft = semilogx(v_frq_stft, mag2db(abs(fftR(v_ir_stft)))  -3*y_sep, mr_stft  , 'Color', cl_stft,  'LineWidth', lw_sm); 
h_sm   = semilogx(v_frq_stft, mag2db(abs(fftR(v_ir_sm)))    -4*y_sep, mr_sm    , 'Color', cl_sm,    'LineWidth', lw_sm);


xlim(xlim_pl)
ylim([-30 10])
set(gca, 'Ytick', [-25:5:5],...
         'YTickLabel', { '' '-20' '' '-10' ''  '0' ''},...
         'XTick', [0.1 0.4 1 2 4 6 8 12 16] * 1000,...
         'XTickLabel', {'0.1', '0.4', '1', '2', '4', '6', '8', '12', '16'});
grid on
box on
legend([h_an h_semi h_sm h_bt h_dft h_stft], ...
       {'Anech.', 'Semi-Anech.', 'Smooth.', 'BT', 'FDT-DFT', 'FDT-STFT'},...
       'Location', 'EastOutside')
ylabel('|H| [dB]')
title('Processed Amplitude Responses')

subplot(312)

y_sep = 2; % Offset between lines, for better display
semilogx(v_frq_stft, mag2db(abs(fftR(v_ir_semian))) - mag2db(abs(fftR(v_ir_anec))), ...
    mr_semi    , 'Color', cl_semi,    'LineWidth',   lw_sm, 'Visible', 'off'); hold on

semilogx(v_frq_stft, mag2db(abs(fftR(v_ir_bt )))    - mag2db(abs(fftR(v_ir_anec))) + 1*y_sep,...
    mr_bt    , 'Color', cl_bt,    'LineWidth', lw_sm);
semilogx(v_frq_stft, mag2db(abs(fftR(v_ir_dft)))    - mag2db(abs(fftR(v_ir_anec))) + 0*y_sep, ...
    mr_dft   , 'Color', cl_dft,   'LineWidth', lw_sm); hold on
semilogx(v_frq_stft, mag2db(abs(fftR(v_ir_stft)))   - mag2db(abs(fftR(v_ir_anec))) - 1*y_sep, ...
    mr_stft  , 'Color', cl_stft,  'LineWidth', lw_sm); hold on
semilogx(v_frq_stft, mag2db(abs(fftR(v_ir_sm  )))   - mag2db(abs(fftR(v_ir_anec))) - 2*y_sep,...
    mr_sm    , 'Color', cl_sm,    'LineWidth', lw_sm);


xlim(xlim_pl)
ylim([-6 4])
set(gca, 'Ytick', [-4:4],...
         'YTickLabel', { '  -4' '' '  -2' '' '   0' ''  '  2' '' ''},...
         'XTick', [0.1 0.4 1 2 4 6 8 12 16] * 1000,...
         'XTickLabel', {'0.1', '0.4', '1', '2', '4', '6', '8', '12', '16'});
grid on
box on
ylabel('\Delta|H| [dB]')
title('Deviation from Anechoic Amplitude')

subplot(313)

y_sep = pi/8;
semilogx(v_frq_stft, unwrap( angle(fftR(v_ir_anec)) -  angle(fftR(v_ir_semian)) ),...
     mr_semi    , 'Color', cl_semi,    'LineWidth',   lw_sm, 'Visible', 'off'); hold on
semilogx(v_frq_stft, unwrap( angle(fftR(v_ir_anec)) -  angle(fftR(v_ir_bt )) ) + 2*pi   + 1*y_sep, ...
    mr_bt    , 'Color', cl_bt,    'LineWidth', lw_sm);
semilogx(v_frq_stft, unwrap( angle(fftR(v_ir_anec)) -  angle(fftR(v_ir_dft)) ) + 2*pi   + 0*y_sep,...
    mr_dft   , 'Color', cl_dft,   'LineWidth', lw_sm); hold on
semilogx(v_frq_stft, unwrap( angle(fftR(v_ir_anec)) -  angle(fftR(v_ir_stft )) ) - 2*pi - 1*y_sep, ...
    mr_stft  , 'Color', cl_stft,  'LineWidth', lw_sm); hold on
semilogx(v_frq_stft, unwrap( angle(fftR(v_ir_anec)) -  angle(fftR(v_ir_sm  )) )         - 2*y_sep, ...
     mr_sm    , 'Color', cl_sm,    'LineWidth', lw_sm);
xlim(xlim_pl)
ylim([-1 1])
set(gca, 'Ytick', [-pi/4 : pi/16 : pi/4],...
         'YTickLabel', { '-pi/4' '' '-\pi/8' '' '0' ...
                         '' '\pi/8' '' '\pi/4'},...
         'XTick', [0.1 0.4 1 2 4 6 8 12 16] * 1000,...
         'XTickLabel', {'0.1', '0.4', '1', '2', '4', '6', '8', '12', '16'});
grid on
box on
ylabel('\Delta\angleH [rad]')
title('Deviation from Anechoic Phase')
pause(.1)

xlabel('f [kHz]')

