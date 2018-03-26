function M_Wsm_MM = csmooth_init_FDTeq( S_cfg )
% M_Wsm_MM = csmooth_init_FDTeq( S_cfg )
%
% initialization of spectral smoothing matrix that is equivalent to
% frequency dependent windowing (FDT)
%
% Modified and specialized implementation of 
%   Hatziantoniou, Panagiotis D., and John N. Mourjopoulos. 
%   "Generalized fractional-octave smoothing of audio and acoustic responses." 
%   Journal of the Audio Engineering Society 48.4 (2000): 259-280.
%
%
% Frequency domain implementation, for half complex spectra!
%
% Input parameters, entries in S_cfg
%   v_frq    : Frequency vector in Hz assiciated to spectra being processed
%              Range: 0 .. fs/2
%   s_win    : Smoothing window shape
%              'rect'           - Rectangular window 
%   s_fscale : < Optional >
%              'FDT_eq': Frequency dependent windowing equivalent
%
%   v_trwin  : truncation window lengths in seconds
% 
% Output Argument:
%   M_Wsm_MM: Smoothing Matrix for use in csmooth_proc, or simple matrix
%   multiplication M_Wsm * v_data, yields complex smoothed spectrum.
%   Adaption of Smoothing Matrix from Eq. (8) for efficient implementation
%   of Eq. (49
%
%
% Author: Florian Denk, January 2017
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


%% Parse input and set defaults
assert(isfield(S_cfg, 'v_frq'),'Missing Input: v_frq')

if ~isfield(S_cfg,'s_fscale')
    S_cfg.s_fscale = 'FDT_eq';
end

if  strcmpi(S_cfg.s_fscale,'FDT_eq') && (length(S_cfg.v_trwin) ~= length(S_cfg.v_frq))
    error('Frequency and Truncation window must have identical lengths')
end
    
    
n_fftR = length(S_cfg.v_frq);
%% Compute Smoothing Index for each frequency

% Compute Smoothing Bandwidth over Frequency
switch S_cfg.s_fscale
    case 'FDT_eq'
        for f = 1:length(S_cfg.v_frq)
            v_Pf = 1./S_cfg.v_trwin; % Equivalent bandwidth in Hz: inverted length in seconds
        end
    otherwise
        error('Invalid Frequency Scaling')
end
    

% Convert to smoothing index, equation (9)
v_m = round(0.5 * v_Pf./ mean(diff(S_cfg.v_frq))); 

%% Compute Windows and fill into Smoothing Matrix
% Implementation of Eq. (49) by Matrix multiplication, much more efficient
% in Matlab than Picking from Smoothing Matrix in Loop
% Also: Generalization to different window classes, not only tapered window
% - substituton of Eq. (5)

M_Wsm_MM = zeros(n_fftR);
for k = 1:n_fftR
    % Window calculation: separate function
    M_Wsm_MM(k,:) = calc_win(n_fftR,...
                             k,...
                             v_m(k),...
                             S_cfg.s_win);
end

end

%% Window calculation
function v_win = calc_win(N,k,m,s_win)
% Compute window as row vector. Window normalized to 1 (sum over entries)
%   N : window length (= number of frequency bins in spectrum)
%   k : window location (= current frequency bin)
%   m : smoothing parameter (= window size)
%   s_win: window type
% 
% Florian Denk, January 2017


v_win = zeros(1,N);

N_winS = 2*m + 1; % Width of desired smoothing range

% Compute Smoothing range indices in frequency vector
v_krange = (k-m) : (k+m);
% Cannot k - v_mrange cannot be smaller than 1 or greater than N
v_krange = v_krange(v_krange > 0);
v_krange = v_krange(v_krange <= N );

% Compute range in temp window
v_wrange = 1:(N_winS);
% Truncate original window in the same manner as original mrange
v_wrange = v_wrange(ismember((k-m) : (k+m) , v_krange) ); 

switch s_win
    case 'rect'
        % Rectangluar window
        v_win( v_krange ) = ones(1,length(v_wrange)) / N_winS;
                
    otherwise
        error('invalid smoothing window shape')
end

% Normalize
v_win = v_win ./ sum(v_win);
end
