function plot_FDTsummary(v_t_ir, M_IR, S_FDwinSteps, n_fft, C_legend, v_flim, ...
    v_tlim, v_fticks, v_spec_lim, plot_irTr, plot_logT, M_colors)
% plot_FDTsummary(v_t_ir, M_IR, S_FDwinSteps, n_fft, C_legend, v_flim, ...
%     v_tlim, v_fticks, v_spec_lim, plot_irTr, plot_logT, M_colors)
%
% Displays summary plot of STFT-based FDT processing.
% Input arguments nr. 6-12 are optional.
%
% For detailed usage, see EXPI_SIMULATIONS and EXPII_MEASUREDDATA
%
% Florian Denk, November 2017
% Dept. Medical Physics & Acoustics, CvO Uni Oldenburg

%% This code is supplementary material of the publication
%
% Denk, F., Kollmeier, B and Ewert, S.D.
% "Removing reflections in semianechoic impulse responses by 
%  frequency-dependent truncation "
% Journal of the Audio Engineering Society 66(3), p. 146-153, 2018
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


%% Parse Inputs
black =       [0 0 0];
grey1 = .45 .* [1 1 1];
grey2 = .75 .* [1 1 1];

if nargin < 6
    v_flim = [0 8000];
end
if nargin < 7
    v_tlim = [-10 60]/1000;
end
if nargin < 8
    v_fticks = [50 200 400 1000 2000 4000];
end
if nargin < 9
    v_spec_lim = [-20 5];
end
if nargin < 10
    plot_irTr = 1;
end
if nargin < 11
    plot_logT = 0;
end
if nargin < 12
    % Define Colors to Plot
    M_colors = [0 0 0;... black, Original
                0 0 1;... blue, FDT
                1 0 0]; % red    TDT
end


% Nearest neighbour interpolation to logarithmic axis: No not change resolution 
s_interp_method = 'nearest';

% Write frequency axis - also logarithmic one, for depiction later
n_frq = n_fft/2+1;
v_frq = linspace(0, S_FDwinSteps.srate/2, n_frq);

v_frq_log = logspace(log10(v_frq(2)/2),log10(v_flim(2)), n_frq);

% Determine bins used for ticks later
for k = 1:length(v_fticks)
    [~,v_inds_fticks(k)] = min( abs(v_frq_log - v_fticks(k)) );
end

% Determine bins for frequency axis
v_inds_flim(1) = find( v_frq_log > v_flim(1), 1, 'first' ) ;
v_inds_flim(2) = find( v_frq_log < v_flim(2), 1, 'last'  ) ;

% F-Ticklabels
c_ftcklabel =  vec2cell( v_fticks );
for ii = 1:length(c_ftcklabel)
    s_tmp = c_ftcklabel{ii};
    if length(s_tmp)>3
    if strcmpi(s_tmp(end-2:end),'000')
        s_tmp = [s_tmp(1:end-3) 'k'];
    end
    end
    c_ftcklabel{ii} = s_tmp;
end


% Compute truncation length of of TD-windowed IR
% Is always minimum of FDT times
v_TD_trncte = min(S_FDwinSteps.v_trwin_length)* ones(n_frq,1);



% Figure Sizes and Ax widths: 

Fontsize = 9;

ax_base = 2.3;

margin_l_cm = 1.2;
margin_u_cm = 1.2;

fig_size = [8.2 margin_u_cm + 3*ax_base + double(plot_logT)*ax_base + 0.2];

margin_l = margin_l_cm /fig_size(1);
margin_u = margin_u_cm /fig_size(2);

plT_w = 2*ax_base  / fig_size(1);% Width of TD plot
plT_h =   ax_base  / fig_size(2);% Height of Time plot

plS_w = 2*ax_base  / fig_size(1);% Width STFT-plot
plS_h = 2*ax_base  / fig_size(2); % Height of STFT plot

plF_w =   ax_base  / fig_size(1);% Width of FD plot
plF_h = 2*ax_base  / fig_size(2);% Height of FD plot



% Interpolate STFT-matrix to logarithmic axis
M_spec = interp1(S_FDwinSteps.v_frq_STFT, mag2db(abs(S_FDwinSteps.M_STFT_orig)), v_frq_log , s_interp_method);
M_spec = M_spec - max(M_spec(:)); % Bring maximum to 0

% Also interpolate truncation length
v_FDwinlength = interp1(S_FDwinSteps.v_frq_STFT, S_FDwinSteps.v_length_trwin_eff, v_frq_log , s_interp_method);

% Calculate spectra
M_Resp = mag2db(abs(fftR(M_IR,n_fft)));
% Interpolate to same logarithmic frequency axis
M_Resp = interp1(v_frq,  M_Resp, v_frq_log );


% Fixed number of integer ticks in Frequency response plot
n_ticks = 7;
MaxTick = floor((1000*v_tlim(2)-0.2) / (n_ticks-1)) * (n_ticks-1);
v_xtick = linspace(0, MaxTick , n_ticks);

% Get FDT windows
M_windows_TD = FDT_windows(S_FDwinSteps);

% Generate figure
h_f = figure('Units','centimeters','position',[10 10 fig_size]);

% plot 1: STFT
h_axS = axes;
% Interpolate STFT Data to logarithmic values of frequency axis -
% logarithmic values!
% % Shift x-axis, so that truncation is displayed to edges of bins
imagesc(S_FDwinSteps.v_t_BlockCentreTimes + S_FDwinSteps.n_shift/2/S_FDwinSteps.srate , ...
    1:n_frq, M_spec, [-100, 0]); hold on

% Colormap: Bone
v_cmap = flipud(bone);
colormap(v_cmap(3:end-4,:));

% plot truncation lines - always
plot(v_FDwinlength , 1:n_frq, '-', 'LineWidth', 2, 'Color', grey1) ;
plot(v_TD_trncte , 1:n_frq, '-.', 'LineWidth', 2, 'Color', grey2);

ylabel('f [Hz]','FontSize',  Fontsize)
xlabel('t [ms]','FontSize' ,  Fontsize)

% Style axis
set(h_axS,'position', [margin_l ...
                       margin_u...
                       plS_w ...
                       plS_h],...
          'FontSize',  Fontsize,...
          'Box', 'on',...
          'XGrid', 'on',...
          'YGrid','On',...
          'ydir','normal', ...
          'ylim',v_inds_flim,...
          'YTick', v_inds_fticks,...
          'YTickLabel', c_ftcklabel,...
          'Xlim', v_tlim,...
          'Xtick', v_xtick/1000,...
          'XTickLabel', vec2cell(v_xtick));
      
% plot 2: IR
h_axT = axes(); % Style later, due to repositioning stuff after legends etc.
h_ir(1) = plot(h_axT, v_t_ir, M_IR(:,1), 'LineWidth', 1.5, 'Color', M_colors(1,:)); hold on % Plot Original
if size(M_IR,2) >= 2
    h_ir(2) = plot(h_axT, v_t_ir, M_IR(:,2), '-','LineWidth', .9, 'Color', M_colors(2,:));
    % Do not plot third input here
end


% Style Axis
set(h_axT,'position', [margin_l ...
                       margin_u + plS_h ...
                       plT_w ...
                       plT_h],...
          'FontSize', Fontsize,...
          'Box', 'on',...
          'XGrid', 'on',...
          'YGrid','On',...
          'XTickLabel','',...
          'XTick', v_xtick/1000,...
          'Xlim', v_tlim,...
          'Ylim', [-round(max(abs(M_IR(:))),1)-.1 round(max(abs(M_IR(:))),1)+.1]);
ylabel('h [a.u.]','FontSize', Fontsize)
set(gca, 'Ytick',[-.6 :.2:.6],...
       'YtickLabel', [{''},vec2cell([-.4 :.2:.4]) ,{''}]);
     
% Also plot effective TD windows
if plot_irTr
v_lim = get(h_axT, 'Ylim');
h_wins = plot(v_t_ir, M_windows_TD(:,1:length(v_t_ir)) .* ( .9* v_lim(2) ) , ...
    'LineWIdth', 1, 'Color', grey1); clear v_lim
T_bt = mean (v_TD_trncte);
n_irwin = 24;
i_end = find(v_t_ir<=T_bt, 1, 'last');
v_btwin = [ones(i_end-n_irwin,1); win_out(n_irwin);zeros(length(v_t_ir)-i_end, 1)];
h_btwin = plot(v_t_ir, v_btwin .* 0.9*max(get(gca, 'Ylim')), '-.', 'LineWidth', 1.5, 'Color', grey2);
uistack(h_wins, 'bottom')
uistack(h_btwin, 'bottom');
end



% plot 3: Frequency Response:
h_axF = axes;
h_ir(1)=  plot(h_axF, M_Resp(:,1), 1:n_frq, 'LineWidth', 1.5, 'Color', M_colors(1,:)); hold on
if size(M_IR,2) >= 2
    h_ir(2)= plot(h_axF, M_Resp(:,2),  1:n_frq, '-','LineWidth', 0.9, 'Color', M_colors(2,:));
end
if size(M_IR,2) >= 3
    h_ir(3)= plot(h_axF, M_Resp(:,3),  1:n_frq, ':', 'LineWidth', 1.5, 'Color', M_colors(3,:));
end
% Put in legend
[h_leg] = legend(h_ir, C_legend, ...
             'Location', 'SouthEastOutside');
% Style axis
set(h_axF,'position', [margin_l + plS_w ...
                       margin_u...
                       plF_w ...
                       plF_h],...
          'FontSize',  Fontsize,...
          'Box', 'on',...
          'XGrid', 'on',...
          'YGrid','On',...
          'XDir', 'reverse', ...
          'YLim',v_inds_flim,...
          'YTick', v_inds_fticks,...
          'YTickLabel', '');
if nargin >= 9
    v_tck_tmp = [-100:5:100];
    c_tck_tmp = vec2cell(v_tck_tmp);
    for cc = 2:2:length(c_tck_tmp)
        c_tck_tmp{cc} = '';
    end
    set(gca, 'XLim', v_spec_lim,...
             'XTick',v_tck_tmp,...
             'XTickLabel', c_tck_tmp);
end
xlabel('|H| [dB]','FontSize',  Fontsize)

% plot 4: log-IR
if plot_logT
h_axT = axes();
h_ir(1) = plot(h_axT, v_t_ir, mag2db(abs(M_IR(:,1))), 'LineWidth', 1.5, 'Color', M_colors(1,:)); hold on % Plot Original
if size(M_IR,2) >= 2
    h_ir(2) = plot(h_axT, v_t_ir, mag2db(abs(M_IR(:,2))), '-','LineWidth', 0.9, 'Color', M_colors(2,:));
end
if size(M_IR,2) >= 3
    h_ir(3) = plot(h_axT, v_t_ir, mag2db(abs(M_IR(:,3))), '--','LineWidth', 0.9, 'Color', M_colors(3,:));
end
% Style Axis
set(h_axT,'position', [margin_l ...
                       margin_u + plS_h + plT_h ...
                       plT_w ...
                       plT_h],...
          'FontSize', Fontsize,...
          'Box', 'on',...
          'XGrid', 'on',...
          'YGrid','On',...
          'XTick', v_xtick/1000,...
          'XtickLabel', '',...
          'Xlim', v_tlim,...
          'Ylim', [-80 0]);
ylabel('log(|h|) [dB]','FontSize', Fontsize)
set(gca, 'Ytick',[-70:10:-10]);
end

% Put in top right corner of STFT
pos_leg = h_leg.Position;
pos_axS = h_axS.Position;
h_leg.Position = [pos_axS(1) + pos_axS(3) - pos_leg(3),...
                  pos_axS(2) + pos_axS(4) - pos_leg(4),...
                  pos_leg(3:4)];

% Include Truncation Lengths
annotation('textbox',[margin_l + plS_w + .05/fig_size(1),...
                      margin_u + plS_h + .05/fig_size(2),...
                      plF_w - .05/fig_size(1),...
                      plF_h - .05/fig_size(2)],...
           'String',{ 'Truncation lengths',...
                      '[ Hz   :  ms]' ,...
                     ['0-133 : '  num2str(round(S_FDwinSteps.v_length_trwin_eff(1)*1000,1))], ...
                     ['  -400 : ' num2str(round(S_FDwinSteps.v_length_trwin_eff(2)*1000,1))],...
                     ['  -667 : ' num2str(round(S_FDwinSteps.v_length_trwin_eff(3)*1000,1))], ...
                     ['  -933 : ' num2str(round(S_FDwinSteps.v_length_trwin_eff(4)*1000,1))], ...
                     [' >933 : '  num2str(round(S_FDwinSteps.v_length_trwin_eff(5)*1000,1))]},...
           'FontSize', Fontsize, 'FontName', 'Times');
end

%%% Auxiliary Functions %%%
function [ M_windows ] = FDT_windows( S_Steps )

N_win  = size(S_Steps.M_WindowMask, 2); % number of windows
N_bins = size(S_Steps.M_WindowMask, 1); % Number of frequency bins


M_windows = zeros(N_bins, round(N_win*S_Steps.n_shift + S_Steps.n_block));

% Analysis-resynthesis window
v_win_ansyn =   S_Steps.v_win_analysis.^2 .* (S_Steps.n_shift / sum(S_Steps.v_win_analysis.^2));

% Filling up

for f = 1:N_bins
    idx = 1;
for n = 1:N_win
    % Fill them - column-wise
    M_windows( f, idx:idx+S_Steps.n_block-1) = M_windows( f, idx:idx+S_Steps.n_block-1) + ...
        S_Steps.M_WindowMask(f, n) .* v_win_ansyn';
    idx = idx + S_Steps.n_shift;
end
end

% Remove zeros from padding
M_windows = M_windows(:, 1 + S_Steps.n_zeros_begin : end);

end

function c_out = vec2cell( v_in)
% c_out = vec2cell( v_in)
%
% Converts Vector to String cell, e.g. for Legend usage


c_out = strread(num2str(v_in), '%s')';
end

