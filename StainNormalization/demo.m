%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Demonstration of a variety of stain normalization methods.
%
% Adnan Khan 
% Department of Computer Science, 
% University of Warwick, UK.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Clear all previous data
clc, clear all, close all;
%mex colour_deconvolution.c;

%% Display results of each method
verbose = 1;

%% Load input & reference image
Source = 'Source_small.png';
img_src=imread(Source);
ref=imread('Ref.png');

%% Stain Normalization using RGB Histogram Specification Method
disp('Stain Normalization using RGB Histogram Specification Method');
[ NormHS ] = RGBHistSpec( img_src, ref, verbose );
% disp('Press key to continue'); pause;

%% Stain Normalization using Reinhard Method
disp('Stain Normalization using Reinhard Method');
[ NormRH ] = Reinhard( img_src, ref, verbose );
% disp('Press key to continue'); pause;

%% Color Deconvolution using online available code with Standard Stain 
%  Matrix, Original credit for C code to G.Landini
disp('Color Deconvolution using Lindini''s C code');
%[s1, s2, s3] = colour_deconvolution(img_src, 'H&E');
% Show deconvolved channels
%[ Hc, Ec, Bgc ] = Channels2RGB( img_src, s1, s2, s3, [], verbose );
% disp('Press key to continue'); pause;

%% Color Deconvolution using Our Implementation with Standard Stain Matrix
disp(['Color Deconvolution using Our Implementation with Standard Stain'...
    ' Matrix ']);
[ ~, H, E, Bg, ~ ] = Deconvolve( img_src, [], verbose );
% disp('Press key to continue'); pause;

%% Display comparative results of two deconvolution implementations
figure,
%subplot(231); imshow(Hc);   title('H (C code)');
%subplot(232); imshow(Ec);   title('E (C code)');
%subplot(233); imshow(Bgc);  title('Bg (C code)');
subplot(234); imshow(H);    title('H (Our Implementation)');
subplot(235); imshow(E);    title('E (Our Implementation)');
subplot(236); imshow(Bg);   title('Bg (Our Implementation)');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

%% Stain Separation using Image specific Stain Matrix for H&E 
disp(['Color Deconvolution using Our Stain matrix estimated using the '...
    'proposed method']);
stain_matrix = EstStainUsingMacenko( img_src );
disp(['Color Deconvolution using Our Implementation with Standard ', ...
    'Stain Matrix']);
[ DCh, H, E, Bg ] = Deconvolve( img_src, stain_matrix', verbose );
% disp('Press key to continue'); pause;

%% Stain Normalization using Macenko Method
disp('Stain Normalization using Macenko et al Method');
[NormMM ] = Macenko(img_src, ref, 255, 0.15, 1, verbose);
% disp('Press key to continue'); pause;

%% Stain Normalization using the Proposed Method
disp('Stain Normalization using the Proposed Method');

% if exist(['normalised/', Source], 'file') ==2
%     delete(['normalised/', Source]);
% end
if exist('normalised/', 'dir') == 7
    rmdir('normalised', 's');
end

if exist('normalised/', 'dir')==0
    mkdir('normalised/');
end

dos(['ColourNormalisation.exe BimodalDeconvRVM filename.txt', ...
    ' Ref.png HE.colourmodel']);
% pause(4);
NormDM = imread(['normalised/', Source]);

if verbose
    figure,
    subplot(131); imshow(ref);   title('Reference Image');
    subplot(132); imshow(img_src);   title('Source Image');
    subplot(133); imshow(NormDM);     title('Normalized (Proposed)');
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
end

% if exist(['normalised/', Source], 'file') ==2
%     delete(['normalised\', Source]);
% end
pause(1);
if exist('normalised/', 'dir') == 7
    rmdir('normalised', 's');
end
% disp('Press key to continue'); pause;

%%  Comparitive Results
disp(' Now Displaying all Results for comparison');
figure,
subplot(231); imshow(ref);          title('Reference');
subplot(232); imshow(img_src);      title('Source');
subplot(233); imshow(NormHS);       title('HistSpec');
subplot(234); imshow(NormRH);       title('Reinhard');
subplot(235); imshow(NormMM);       title('Macenko');
subplot(236); imshow(NormDM);       title('Proposed');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

%%
disp('End of Demo');
