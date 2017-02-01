function [ stain_matrix ] = EstStainUsingAdnan( I, Io, beta, alpha )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deconvolve: Deconvolution of an RGB image into its constituent stain
% channels
%
%
% Input:
% I                 - RGB input image;
% beta              - OD threshold for transparent pixels
% alpha             - tolerance for the pseudo-min and pseudo-max
% Output:
% stain_matrix      - hematoxylin image;
%
% References:
% A method for normalizing histology slides for quantitative analysis. M.
% Macenko et al., ISBI 2009
%
%
% Copyright (c) 2013, Adnan Khan
% Department of Computer Science,
% University of Warwick, UK.
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run in DEMO Mode
if nargin<1
    I = imread('hestain.png');
end



% OD threshold for transparent pixels
if ~exist('beta', 'var') || isempty(beta)
    beta = 0.15;
end

% tolerance for the pseudo-min and pseudo-max
if ~exist('alpha', 'var') || isempty(alpha)
    alpha = 1;
end

% transmitted light intensity
if ~exist('Io', 'var') || isempty(Io)
    Io = 255;
end

if size(I,3)<3, errordlg('Image must be RGB'); return; end

I = reshape(double(I), [], 3);

% calculate optical density
OD = -log((I+1)/Io);

% remove transparent pixels
ODbar = OD(~any(OD < beta, 2), :);

% calculate eigenvectors
[V, ~] = eig(cov(ODbar));

% project on the plane spanned by the eigenvectors corresponding to the two
% largest eigenvalues
THETA = ODbar*V(:,2:3);

PHI = atan2(THETA(:,2), THETA(:,1));

% find the robust extremees (min and max angles) 
minPhi = prctile(PHI, alpha);
maxPhi = prctile(PHI, 100-alpha);

% Bring the extreme angles back to OD Space
VEC1 = V(:,2:3)*[cos(minPhi); sin(minPhi)];
VEC2 = V(:,2:3)*[cos(maxPhi); sin(maxPhi)];

% Make sure that Hematoxylin is first and Eosin is second vector

if VEC1(1) > VEC2(1)
    stain_matrix = [VEC1 VEC2];
else
    stain_matrix = [VEC2 VEC1];
end

% stain_matrix = stain_matrix;
end

