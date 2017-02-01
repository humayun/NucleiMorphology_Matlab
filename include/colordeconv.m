function [iA, iS, iAS, iR, dA, dS, dAS, dR] = colordeconv(I, nstains, varargin)
% This function performs the color deconvolution on the input multi-channel
% image I into NSTAINS stain channels by using the non-negative matrix
% factorization.
%
% Input:
%       I: the input multi-channel image
%       NSTAINS: number of stains
%       VARARGIN: optional parameters
%
% Output:
%       IA: (H * W * NSTAINS) matrix representing the amount of stains in
%       intensity.
%       IS: (NSTAINS * NCHANNELS) matrix representing the absorbance factor of each
%       channel for each stain in intensity.
%       IAS: (NSTAINS * 1) cell, each element of which is a matrix reprensenting the
%       multi-channel image of the corresponding stain.
%       IR: matrix with the same size of the input image, representing the residual
%       intensity image of the color deconvolution.
%       DA: (H * W * NSTAINS) matrix representing the amount of stains in density.
%       DS: (NSTAINS * NCHANNELS) matrix representing the absorbance factor of each
%       channel for each stain in density.
%       DAS: (NSTAINS * 1) cell, each element of which is a matrix reprensenting the
%       multi-channel image of the corresponding stain.
%       DR: matrix with the same size of the input image, representing the residual
%       density image of the color deconvolution.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING, USING OR
% MODIFYING.
%
% By downloading, copying, installing, using or modifying this
% software you agree to this license.  If you do not agree to this
% license, do not download, install, copy, use or modifying this
% software.
%
% Copyright (C) 2010-2010 Baochuan Pang <babypbc@gmail.com>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
% General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% check arguments
    % input image should be color image
    ndimsI = ndims(I);
    if ndimsI ~= 3,
        error('colordeconv:InputImageDimensionError', ...
              'Input image I should be 3-dimensional!');
    end

    % number of stains should not exceed the number of channels
    nchannels = size(I, 3);
    if nstains > nchannels,
        error('colordeconv:InputStainsError', ...
              'Number of stains should not exceed the number of channels!');      
    end
   
%% compute optical density from the intensity image
    % height and width
    h = size(I, 1);
    w = size(I, 2);

    % reshape each pixel as a row of values
    vI = reshape(I, h*w, nchannels);
   
    % compute optical density
    %vD = compute_optical_density(vI);
    
    npixels = size(vI, 1);
    nchannels = size(vI, 2);
    maxIs = max(vI, [], 1);
    logMaxIs = [5.5413, 5.5413, 5.5413]; %log(maxIs);
    logIs = [5.5413, 5.5413, 5.5413] %log(vI);
    vD = ones(npixels, nchannels) * diag(logMaxIs) - logIs;
   
%% compute the amount of stain and absorption factor of stains
    % non-negative matrix factorization, where As denotes the
    % amount of stains in each pixel, and S indicates the
    % absorption factor of each stain
    [vdA, dS] = nnmf(vD, nstains);
   
%% reshape to the original image size and reconstruct to
    %  intensity
    dA = reshape(vdA, h, w, nstains);   % density of A
    iA = exp(-dA);              % intensity of A
    iS = exp(-dS);              % intensity of S
   
%% calculate the color image for each stain
    dAS = cell(nstains, 1);
    iAS = cell(nstains, 1);
    for i = 1 : nstains,
            vdAS = vdA(:, i) * dS(i, :);
            dAS{i} = reshape(vdAS, h, w, nchannels);
            iAS{i} = exp(-dAS{i});
    end
       
%% calculate the residual
    vdR = vD - (vdA * dS);
    dR = reshape(vdR, h, w, nchannels);
    iR = exp(-dR);
end

