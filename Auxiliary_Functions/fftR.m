function vf_sig = fftR(v_sig,N)
% vf_sig = fftR(v_sig, N)
% input: v_sig: input signal to be transformed, column vector or Matrix
%        N    : length of fft
% Function that does FFT for real signals and truncates spectrum to half
% sanpling frequency. Otherwise same behaviour as matlab's fft function.
% Input can be vector or matrix, for matrix fft will be performed
% column-wise.
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



if nargin < 2
    N = size(v_sig,1);
end
if ~isreal(v_sig)
    warning('This function should only be applied to real-valued signals')
end

% do FFT
vf_sig = fft(v_sig,N);

% truncate spectrum
vf_sig = vf_sig( 1 : ceil(end/2)+1,:,:,:);

end