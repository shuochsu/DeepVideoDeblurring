% TV-L1 Optical Flow Estimation--
% Javier Sánchez Pérez, Enric Meinhardt-Llopis, Gabriele Facciolo
%
% @article{ipol.2013.26,
%    title   = {{TV-L1 Optical Flow Estimation}},
%    author  = {Sánchez Pérez, Javier and Meinhardt-Llopis, Enric and Facciolo, Gabriele},
%    journal = {{Image Processing On Line}},
%    volume  = {3},
%    pages   = {137--150},
%    year    = {2013},
%    doi     = {10.5201/ipol.2013.26},
% }
%
%
% MATLAB MEX function for tvl1optflow estimation. The code is provided by
% Javier Sanchez Peres <jsanchez@dis.ulpgc.es>
% The MEX file was created by Mauricio Delbracio <mdelbra@gmail.com> and
% the optical flow estimation was modified to be computed in a 1/3 smaller
% image and then it is upsampled to the original resolution
% 
%
% This program is free software: you can use, modify and/or redistribute it
% under the terms of the simplified BSD License. You should have received a
% copy of this license along this program. If not, see
% <http://www.opensource.org/licenses/bsd-license.html>.
%
% Copyright (C) 2011, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
% All rights reserved.
%
%
% USAGE:
%
% flow = tvl1flow(I0,I1,zfini,tau,lambda,theta,nscales,zfactor,nwarps,epsilon,verbose);
%
% Example parameters (more info in the IPOL website) 
%  zfini=0.3 (Initial subsampling)
%  tau=0.25;
%  lambda = 0.1; (Most important parameter, regularization)
%  theta=0.3;
%  nscales=4;
%  zfactor=0.5;
%  nwarps=5;
%  epsilon=0.01
%  verbose=1;
%
