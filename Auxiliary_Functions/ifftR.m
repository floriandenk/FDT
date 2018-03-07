function v_sig = ifftR(vf_sig)
% v_sig = ifftR(vf_sig)
% Function that does inverse FFT for truncated spectra of real-valued
% time-domain signals
% Input can be vector or matrix, for matrix fft will be performed
% column-wise (same behaviour as matlab's fft function)
%
% Florian Denk, March 2016
% Dept. Medical Physics & Acoustics, CvO Uni Oldenburg

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


% Fill other half of spectrum with complex-conjugate
if mod(size(vf_sig,1), 2)
    vf_sig = [vf_sig;conj(vf_sig(end-1:-1:2,:,:,:))];
else
    vf_sig = [vf_sig;conj(vf_sig(end-2:-1:2,:,:,:))];
end

% Perform inverse transform, and discard imaginary parts
v_sig = real(ifft(vf_sig));
end