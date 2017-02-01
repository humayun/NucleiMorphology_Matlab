
function [Norm ] = Macenko(Source, Target, Io, beta, alpha, verbose)
%
%
% Macenkono: Normalize the appearance of an RGB Source Image to the
% Reference RGB image.
% 
%
% Input:
% Source    - RGB Source image;
% Target    - RGB Reference image;
% A        - (optional) transmitted light intensity (default 255)
% beta      - (optional) OD threshold for transparent pixels (default 0.15)
% alpha     - (optional) tolerance for the pseudo-min and pseudo-max
%                        (default 1)
% verbose   - (optional) Display Results or not?
%             Default value = 0, don't display 
%
%
% Output:
% Norm      - Normalized RGB image
%
% References:
% A method for normalizing histology slides for quantitative analysis. M.
% Macenko et al., ISBI 2009
%
%
% Adnan Khan
% Department of Computer Science,
% University of Warwick, UK.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


% Run in DEMO Mode
if nargin<1
    Source = imread('hestain.png');
    Target = imread('hestain.png');
end

if ~exist('verbose', 'var') || isempty(verbose)
   verbose = 0; 
end

% transmitted light intensity
if ~exist('A', 'var') || isempty(A)
    A = 255;
end

% OD threshold for transparent pixels
if ~exist('beta', 'var') || isempty(beta)
    beta = 0.15;
end

% tolerance for the pseudo-min and pseudo-max
if ~exist('alpha', 'var') || isempty(alpha)
    alpha = 1;
end

%%
[h, w, ~] = size(Source);

% Estimate Stain Matrix for Taret Image
Mtarget = EstStainUsingMacenko(Target, A, beta, alpha);

% Perform Color Deconvolution of Target Image to get stain concentration
% matrix
[ C,~,~,~,Mtarget ] = Deconvolve( Target, Mtarget' );

% Vectorize to 3 x N matrix
C = reshape(C, [], 3)';

% Find the 99th percentile of stain concentration (for each channel)
maxCTarget = prctile(C, 99, 2);    


%% Repeat the same process for input/source image

% Estimate Stain Matrix for Source Image
MSource = EstStainUsingMacenko(Source, A, beta, alpha);

% Perform Color Deconvolution of Source Image to get stain concentration
% matrix
[ C,~,~,~,~ ] = Deconvolve( Source, MSource' );

% Vectorize to 3 x N matrix
C = reshape(C, [], 3)';

% Find the 99th percentile of stain concentration (for each channel)
maxCSource = prctile(C, 99, 2);


%% MAIN NORMALIZATION STUFF
% 
C = bsxfun(@rdivide, C, maxCSource);
C = bsxfun(@times,   C, maxCTarget);


%% VISUALIZATION
% Reconstruct the RGB image 
Norm = A*exp(-Mtarget * C)';
Norm = reshape(Norm, h, w, 3);
Norm = uint8(Norm);

% Display results if verbose mode is true
if verbose==1;
    figure; 
    subplot(131); imshow(Target);   title('Reference Image');
    subplot(132); imshow(Source);   title('Source Image');
    subplot(133); imshow(Norm);     title('Normalized (Macenko)');
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
end

end


