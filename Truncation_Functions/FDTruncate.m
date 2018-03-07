function [v_ir_proc, S_intermediate] = ...
                    FDTruncate(v_ir_in, ...
                               v_trwin_length, ...
                               srate, ...
                               n_block, ...
                               n_shift,...
                               v_peakrange)
% function [v_ir_proc, S_intermediate] = ...
%                     FDTruncate(v_ir_in, ...
%                                v_trwin_length, ...
%                                srate, ...
%                                n_block, ...
%                                n_shift, ...
%								 v_peakrange)
%
% Performs Frequency-Dependent Truncation of acoustic impulse response
% over a Short Time Fourier transform
% 
% For details on the utilization, see EXPI_SIMULATIONS
%
%
% Inputs: 
%   v_ir_in:        Impulse Response to be truncated, as column vector
%   v_trwin_length: Truncation lengths for each frequency bin in seconds. 
%                   Length must be n_block/2+1
%   srate:          Sampling Rate of IR
%   n_block:        Blocklength of STFT analysis, i n samples. Recommended: ca. 3 ms,
%   n_shift:        Hop Size in STFT analysis in samples. n_block/2 is
%                   sufficient. Smaller hop sizes result in a smoother
%                   Spectrogram representation, however the effective fade
%                   out length in each bin is increased (on the other hand,
%                   this does not do much harm)
%   v_peakrange:    (Optional)
%                   Index range in v_ir_in, where the peak of the IR is 
%                   searched. This should be specified in cases where the
%                   amplitude of a reflection may be higher than the
%                   initial peak.
%                   
% Outputs:
%   v_ir_proc:      Processed impulse response
%   S_intermediate: Intermediate variables for plotting etc.:
%      M_STFT_orig:           STFT matrix prior to truncation(Frequency in
%                             1st dimension)
%      v_frq_STFT:            Frequency vector of STFT blocks
%      v_win_analysis:        Analysis (= Reconstruction) window
%      n_block:               Block length of STFT Transform
%      n_shift:               Hop Size of STFT Transform
%      n_zeros_begin:         Zeros appended to beginning of IR prior to STFT
%      n_delay_ir:            Position of IR peak (including padded zeros)
%      srate:                 Sampling rate
%      v_trwin_length:        Desired Truncation window lengths
%      v_length_trwin_eff:    Effective Truncation window lengths, in STFT domain
%      v_t_BlockCentreTimes:  STFT Block centre times wrt IR peak, including
%                             padded zeros
%      v_t_BlockEndTimes:     STFT Block end times wrt IR peak, including
%                             padded zeros
%      M_WindowMask:          Binary Truncation Mask applied to STFT
%      M_STFT_final:          STFT representation after truncation
% 
% IR is converted to the STFT domain (with blocksize n_block and shift 
% n_shift, both in samples), and for each bin a different window
% length (v_trwin_lengh, starting from IR peak) is applied by means of 
% a binary mask. 
%
% Prior to the STFT transform zeros are appended to the
% original IR in the beginning such that one frame is ending exactly at the
% shortest specified truncation length (usually in the high-frequency
% domain). Truncation is performed that the specified duration marks the
% point where the fade-out window in this frequency bin is exactly 0.
% Effective truncation lengths, depending on the STFT block positions, are
% contained in the S_intermediate Struct
%
%
% Florian Denk, March-December 2017
% Dept. Medical Physics and Acoustics, CvO Uni Oldenburg

%% This code is supplementary material of the publication
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



%% Parse input and check incompatibilities

if nargin < 6
    v_peakrange = 1:length(v_ir_in);
end


assert(size(v_ir_in, 2) == 1 && size(v_ir_in, 3) == 1 && size(v_ir_in, 4) == 1,...
       'Input data may only have one dimension')
   
assert(n_block/n_shift >= 2, 'STFT shift may not be larger than half block length')

assert(~mod(n_block, 2), 'Block length shoud be even');

assert(length(v_trwin_length) == n_block/2+1 || isscalar(v_trwin_length), ...
    'Truncation window length vector must match block length or be scalar')

% Check for length of desired windows
t_ir = size(v_ir_in,1)/srate;
if max(v_trwin_length(:)) > t_ir
    warning('In at least one frequency bin the orginial IR is shorter than the truncation length')
    v_trwin_length = min(v_trwin_length, t_ir);
end

%% Analysis
% Append zeros at beginning to avoid boundary effects in STFT processing
v_ir_pd = [zeros(n_block,1) ; v_ir_in];

% Compute delay in IR
[~,n_delay_ir] = max(abs(hilbert(v_ir_pd(v_peakrange))));

% Index where shortest truncation length ends
n_ind_end = n_delay_ir + round( min(v_trwin_length(:))*srate );
n_ind_clastwin = n_ind_end - n_block/2; % desired center index of last block

% Append more zeros at the beginning, such that full STFT analysis has a
% window ending at shortest truncation length
n_extraZeros = n_shift - mod(n_ind_clastwin - n_shift, n_shift);
v_ir_pd = [zeros(n_extraZeros,1) ; v_ir_pd];

% Include appeded zeros into IR delay
n_delay_ir = n_delay_ir + n_extraZeros;

% Perform STFT
[M_STFT, v_frq, v_t_blocks, v_win_analysis] = mediSTFT(v_ir_pd, n_block, n_shift, srate, 'sqrhann');

% Store intermediate variables
if nargout > 1
    S_intermediate.M_STFT_orig     = M_STFT;
    S_intermediate.v_frq_STFT      = v_frq;
    S_intermediate.v_win_analysis  = v_win_analysis;
    S_intermediate.n_block         = n_block;
    S_intermediate.n_shift         = n_shift;
    S_intermediate.n_zeros_begin   = n_block + n_extraZeros;
    S_intermediate.n_extra_zeros   = n_extraZeros;
    S_intermediate.n_delay_ir      = n_delay_ir;
    S_intermediate.srate           = srate;
    S_intermediate.v_trwin_length   = v_trwin_length;
end

%% Frequency Dependent Truncation

% Reset time vector by IR delay
v_t_blocks = v_t_blocks - (n_delay_ir)/srate; % Now: shift time axis by this
% Time vector, Block end points
v_t_BlockEndTimes = v_t_blocks + n_block/2/srate; 

% Store Time vectors
if nargout > 1
    S_intermediate.v_t_BlockCentreTimes = v_t_blocks;
    S_intermediate.v_t_BlockEndTimes    = v_t_BlockEndTimes;
    
end
                    
% Write windowing matrix
M_WindowMask = zeros(size(M_STFT));
v_length_trwin_eff = zeros(length(v_frq), 1);

% Write desired frames to 1, others to zero. 
for f = 1:length(v_frq)
    M_WindowMask( f, v_trwin_length(f) - v_t_BlockEndTimes > - 1/srate ) = 1;
    % Store effective truncation length: End time of last block wrt IR peak
    v_length_trwin_eff(f) =  v_t_BlockEndTimes( find( M_WindowMask(f,:) == 1, 1, 'last' ));
end

% Store variables
if nargout > 1
    S_intermediate.M_WindowMask = M_WindowMask;
    S_intermediate.v_length_trwin_eff = v_length_trwin_eff;
end

% Apply Windowing Matrix = Perform frequency dependent windowing, eq. (6)
M_STFT = M_STFT .* M_WindowMask;

% Truncate STFT representation - discard frames at end where only 0s are contained
n_frame_max = find(sum(M_WindowMask,1) > 0, 1, 'last');
M_STFT = M_STFT(:,1:n_frame_max);

% Store variables
if nargout > 1
    S_intermediate.M_STFT_final = M_STFT;
end


%% Reconstruction
% Do weighted overlap-add
v_ir_proc = mediWOLA(M_STFT, n_shift, v_win_analysis);

% Remove zeros appended in the beginning
v_ir_proc = v_ir_proc(1 + n_block + n_extraZeros : end);

if nargout > 1
    S_intermediate.M_STFT_final = M_STFT;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary Functions  %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M_STFT, v_frq, v_t, v_win_an] = ...
                        mediSTFT(v_in,n_block,n_shift, srate, s_win)
% [M_out, v_frq, v_t, v_win_syn, v_win_an] = ...
%                         mySTFT(v_in,n_block,n_shift, srate, s_win)
% 
% Performs short-time fourier transform with overlapped windows
%
% Inputs: v_in: input vector, data in first dimension. No passing arrays
%         n_block: Length of block in samples
%         n_shift: Block shift in samples
%         srate:   Sampling rate in Hz, Default: 48000
%         s_win:     Window function. Options
%                 Hann, Hanning, Hamming, Blackman, rect, trian
%                 sqrHann: Square-root hann window
%
% Outputs: M_STFT:    Signal in STFT domain, each column containing the half
%                     spectrum of one frame
%          v_frq:     Frequency vector for M_STFT
%          v_t:       Time vector for M_STFT, giving window centres
%          v_win_an:  Used Analysis Window
%
% Florian Denk, 28.02.2017
% Dept. Medical Physics and Acoustics, Uni Oldenburg

% Parse input
assert(size(v_in, 2) == 1, 'Input vector may only have one dimension');

if nargin < 4
    srate = 48000;
end
if nargin < 5
    s_win = 'sqrHann';
end


% compute number of frames
n_Frames = floor((length(v_in) - n_block) / n_shift) + 1;

% create analysis window
switch lower(s_win)
    case {'hann', 'hanning'}
        v_win_an  = hann(n_block, 'periodic');
    case 'sqrhann'
        v_win_an  = sqrt(hann(n_block, 'periodic'));
    case 'hamming'
        v_win_an  = hamming(n_block, 'periodic');
    case 'blackman'
        v_win_an  = blackman(n_block, 'periodic');
    case 'trian'
        v_win_an = triang(n_block);
    case 'rect'
        v_win_an = ones(n_block, 1);
    otherwise
        error(['Unknown analysis window:' s_win])
end

% generate matrix containing block indeces
matIdx = bsxfun(@plus, (1:n_block).', (0:n_Frames - 1) .* n_shift);

% split signal into frames, apply window function and perform fourier
% transformation
M_STFT = fftR(bsxfun(@times, v_in(matIdx), v_win_an));

% Compute time and frequency vectors, if desired for output
if nargout > 1
    T_shft = n_shift/srate;
    v_t = ( 0:T_shft:T_shft*(n_Frames-1) ) + n_block/2/srate;
    v_frq = linspace(0, srate/2, size(M_STFT, 1));
end
end

function [v_out, v_env] = mediWOLA(M_STFT, n_shift, v_win_an)
% v_out = myWOLA(M_STFT, n_shift, v_win)
%
% Weighted Overlapp-Add for STFT-domain signal M_STFT, as created by
% function mySTFT. 
% 
% Inputs: M_STFT:    Signal in STFT domain, each column containing the half
%                    spectrum of one frame
%         n_shift:   Shift between frames in samples
%         v_win_syn: Analysis window used for obtaining STFT
%
% Florian Denk, 28.02.2017
% Dept. Medical Physics and Acoustics, Uni Oldenburg

% Parse input
n_block  = size(v_win_an,  1); % Block length
n_Frames = size(M_STFT, 2);    % Number of frames

% Calculate corresponding synthesis window with the WOLA constraint
% See Handbook of Speech Processing, Chapter 12
v_win_syn = v_win_an .* (n_shift / sum(v_win_an.^2));

% Initialize output
v_out = zeros((n_Frames-1)*n_shift + n_block, 1);

% Do overlap-add
idx = 1;
for i_fr = 1:n_Frames
    v_out(idx:idx+n_block-1) = v_out(idx:idx+n_block-1) + ...
                               ifftR(M_STFT(:,i_fr)) .* v_win_syn;
    idx = idx + n_shift;
end

% If desired, calculate WOLA envelope created by analysis-reconstruction chain
if nargout>1
    v_env = zeros(size(v_out));
    idx = 1;
    for fr = 1:size(M_STFT, 2)
    v_env(idx:idx+n_block-1) = v_env(idx:idx+n_block-1) + (v_win_an .* v_win_syn);
    idx = idx + n_shift;
    end
end
end


