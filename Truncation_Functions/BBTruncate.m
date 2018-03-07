
function v_ir_out = BBTruncate( v_ir_in,n_start,n_len,n_startramplen,n_endramplen )
% v_ir_out = BBTruncate( v_ir_in, n_start, n_len, n_startramplen, n_endramplen )
% 
% Truncation of impulse responses
% returns cut version with length 'n_len' of signal 'v_ir_in', starting from
% sample defined by 'n_start' faded in by a hann ramp of length 
% 'n_startramplen' and faded out by a hanning ramp of length 'n_endramplen'.
% Input can be a more-dimensional array, data must be in first dimension
%
% Author: Florian Denk, September 2016, based on scripts by Stephan Ewert
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


K = size(v_ir_in); K = K(2:end);
n_end=n_len+n_start-1;

% Generate ramps - half hann windows
v_startwin = win_in (n_startramplen);
v_endwin   = win_out(n_endramplen);

% Cut range, and apply ramps
v_ir_out = [v_ir_in(n_start:(n_start+n_startramplen-1),:,:,:,:)      .* repmat(v_startwin,[1 K]);...  
            v_ir_in((n_start+n_startramplen):(n_end-n_endramplen),:,:,:,:)     ;...
            v_ir_in((n_end-n_endramplen+1):n_end,:,:,:,:)      .* repmat(v_endwin,[1 K])];
    
end
