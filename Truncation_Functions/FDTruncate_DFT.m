
function v_ir_out = FDTruncate_DFT(v_ir_in, v_tr_lengths, srate, v_peakrange)
% v_ir_out = FDTruncate_DFT(v_ir_in, v_tr_lengths, srate)
%
% Performs DFT-based Frequency Dependent Truncation
%
% Inputs:  v_ir_in: Original impulse response as column vector
%          v_tr_lengths: truncation lengths in seconds, spaced equally
%              across the sepctrum
%          srate: sampling rate in Hz
%          v_peakrange: index range, where peak of the IR is searched
%
% For details on the utilization, see EXPI_SIMULATIONS
%
% Florian Denk, 27.5.2017
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

if nargin < 4
  v_peakrange = 1:length(v_ir_in);
end

% Windows for windowing in and out
n_in  = 128;  v_win_in  = win_in(n_in);
n_out =  16;  v_win_out = win_out(n_out);

% Get delay of main peak
[~,n_delay_ir] = max(abs(hilbert(v_ir_in(v_peakrange))));

% Adapt fft length to longst truncation length
n_fft = 2^nextpow2(n_delay_ir + round(max(v_tr_lengths(:))*srate));

% Compute frequency vector for according half-spectrum
n_frq = n_fft/2 + 1;
v_frq = linspace(0, srate/2, n_frq)';

% Interpolare truncation lengths to new frequency vector
v_tr_lengths = interp1(linspace(0, srate/2, length(v_tr_lengths)), ...
                       v_tr_lengths,...
                       v_frq, 'nearest');

% Initial truncation to longest truncation length
if n_fft > length(v_ir_in)
    v_ir_in = [v_ir_in; zeros(n_fft - length(v_ir_in), 1)];
end
v_ir_in = v_ir_in(1:n_fft);
v_ir_in(1:n_in) = v_ir_in(1:n_in) .* v_win_in;
v_ir_in(end-n_out+1:end) = v_ir_in(end-n_out+1:end) .* v_win_out;


% Do FDT-DFT (implementation of eq. (3)
vf_out = zeros(n_frq, 1);
% Loop over frequencies
for k = 1:n_frq
    % Index after which region gets truncated
    n_max = n_delay_ir + round(v_tr_lengths(k) * srate);
    
    % compute according window
    v_win = zeros(n_fft, 1);
    v_win(1:n_in) = v_win_in;
    v_win(n_in+1 : n_max-n_out) = 1;
    n_out = round( (n_max - n_delay_ir)/4 );
    v_win(n_max -n_out + 1:n_max) = win_out( n_out );
    
    % Apply truncation window to impulse response
    v_tmp = v_ir_in .* v_win;

    % Compute DFT sum for the current frequency bin 
    for n = 1:n_max
        vf_out(k) = vf_out(k) + v_tmp(n) * exp(- 1i*2*pi*  ((k-1)*(n-1)) / n_fft);
    end
end

% Transform back to time domain, eq. (4)
v_ir_out = ifftR(vf_out);
