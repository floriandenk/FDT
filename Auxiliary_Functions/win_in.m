function v_win_in = win_in(N)
% v_win_in = win_in(N)
%
% Returns asceting half -hann window of N samples
%
% Florian Denk, March 2016
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


v_win_in = hanning(2*N);
v_win_in = v_win_in(1:end/2,:);