function [ DCh, H, E, Bg, M ] = Deconvolve( I, M, verbose )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deconvolve: Deconvolution of an RGB image into its constituent stain
% channels
% 
%
% 
%
% Input:
% I         - RGB input image;
% M        - (optional) Stain matrix. If stain matrix is not provided, the
%             function uses Ruifork and Johnston's method of stain
%             normalization.
% verbose   - (optional) Display Results or not?
%             Default value = 0, don't display 
%
%
% Output:
% DCh       - Deconvolved Channels concatatenated to form a stack. 
%             Each channel is in double format (instead of uint8)
% E1        - hematoxylin image (RGB);
% E2        - eosin image (RGB);
% E3        - Channel-3 image (RGB) 
%
% Reference:
% [1] Ruifrok AC, Johnston DA. Quantification of histochemical staining by
% color deconvolution. Analytical & Quantitative Cytology & Histology 2001;
% 23: 291-299.
%
%
% Example:
%           I = imread('hestain.png');
%           [ DCh, H, E, Bg, M ] = Deconvolve( I, [], 1);
%
%
%
%
% Adnan Khan
% Department of Computer Science,
% University of Warwick, UK.
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run in DEMO Mode
if nargin<1
    I = imread('hestain.png');
end



%% Sanity check

[h, w, c] = size(I);

% Image must be RGB
if c<3, errordlg('Image must be RGB'); return; end


%% Display results or not?
if ~exist('verbose', 'var') || isempty(verbose)
   verbose = 0; 
end


%% Default Color Deconvolution Matrix proposed in Ruifork and Johnston [1]
if ~exist('M', 'var') || isempty(M)
   M = [   0.644211 0.716556 0.266844; 
           0.092789 0.954111 0.283111; 
       ]; 
end


%% Add third Stain vector, if only two stain vectors are provided. 
% This stain vector is obtained as the cross product of first two
% stain vectors. 
if size (M,1) < 3
    M = AddThirdStainVector(M);
end


%% MAIN IMPLEMENTATION OF METHOD

% the intensity of light entering the specimen (see section 2a of [1])
Io = 255;

% Convert to double
J = double(I);

% Vectorize
J = reshape(J, [], 3);

% calculate optical density
OD = -log((J+1)/Io);
Y = reshape(OD, [], 3)';
% determine concentrations of the individual stains
% M is 3 x 3,  Y is 3 x N, C is 3 x N
% C = pinv(M+eps)*Y;
C = M \ Y;

% Stack back deconvolved channels
DCh = reshape(C', h, w, 3);

%% Output and Visualization
% Generate Output Images
if nargout > 1 || verbose
    H = Io*exp(-M(:, 1) * C(1, :));
    H = reshape(H', h, w, 3);
    H = uint8(H);
end

if nargout > 2 || verbose
    E = Io*exp(-M(:, 2) * C(2, :));
    E = reshape(E', h, w, 3);
    E = uint8(E);
end

if nargout > 3 || verbose
    Bg = Io*exp(-M(:, 3) * C(3, :));
    Bg = reshape(Bg', h, w, 3);
    Bg = uint8(Bg);
end

if verbose,
    figure,
    subplot(141); imshow(I); title('Source');
    subplot(144); imshow(Bg); title('Background');
    subplot(142); imshow(H); title('Hamatoxylin');
    subplot(143); imshow(E); title('Eosin');
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
end

