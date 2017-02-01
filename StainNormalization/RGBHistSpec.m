function [ Norm ] = RGBHistSpec( Source, Target, verbose )

% RGBHistNorm: Normalize the RGB histogram of Source Image with respect to
% a target image
% 
%
% Input:
% Source         - RGB image that need to be normalized;
% Target         - RGB REFERENCE image 
% verbose        - (optional)Display Results (including Histograms).
%                   Default value = 0, don't display histograms
%
% Output:
% Norm           - Normalized RGB image
%
% References:
% A. Jain, Fundamentals of digital image processing. Prentice-Hall, 1989.
%
%
% Copyright (c) 2013, Adnan Khan
% Department of Computer Science,
% University of Warwick, UK.
% 
% 

if ~exist('verbose', 'var') || isempty(verbose)
   verbose = 0; 
end

if ~exist('Source', 'var') || isempty(Source)
   errordlg('Need Reference Image');
   return
end

if ~exist('Target', 'var') || isempty(Target)
   errordlg('Need Reference Image');
   return;
end


% Separate Source image’s color channel
SourceR = Source(:,:,1);
SourceG = Source(:,:,2);
SourceB = Source(:,:,3);

%Separate Target/reference image’s color channel

TargetR = Target(:,:,1);
TargetG = Target(:,:,2);
TargetB = Target(:,:,3);

%% 

% Compute Source image histograms
HnSourceR = imhist(SourceR)./numel(SourceR);
HnSourceG = imhist(SourceG)./numel(SourceG);
HnSourceB = imhist(SourceB)./numel(SourceB);

% Compute Target/reference image histograms
HnTargetR = imhist(TargetR)./numel(TargetR);
HnTargetG = imhist(TargetG)./numel(TargetG);
HnTargetB = imhist(TargetB)./numel(TargetB);

% Histogram specification, using Target/reference image's histogram
NormR = histeq(SourceR,HnTargetR);
NormG = histeq(SourceG,HnTargetG);
NormB = histeq(SourceB,HnTargetB);

% Concatenate Channels
Norm = cat(3, NormR, NormG, NormB);

%% Plot histogram & Display Image
if verbose
    %Show Image
    figure;
    subplot(131); imshow(Target,'InitialMagnification','fit'); title('Reference Image');
    subplot(132); imshow(Source,'InitialMagnification','fit'); title('Source Image');
    subplot(133); imshow(Norm,'InitialMagnification','fit'); title('Normalized Image');
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
end
end

