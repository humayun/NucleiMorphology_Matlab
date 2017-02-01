%%
clear all; close all; clc
addpath('H:\CodingStuff\Matlab\stain_normalization_toolbox'); 

%%
path = 'H:\Datasets\TMA_DCIS_Standford\C-TA-239-00.31.HE_40x_4\';
imageExt = '.tiff';
srcFiles = dir(strcat(path,'*',imageExt));
verbose = 0;    % Display results of each method

%% Target Image Stain Estimation
Target = imread('H:\CodingStuff\Matlab\stain_normalization_toolbox\ref2.tiff');

%%
for i = 1:length(srcFiles)
    [~, imageName, imageExt] = fileparts(srcFiles(i).name);
    im = imread( strcat(path, imageName, imageExt)) ;

    disp(['Color Deconvolution with Standard Stain Matrix ']);
    [ ~, H, E, ~, ~ ] = Deconvolve( im, [], verbose );

    disp(['Color Deconvolution using Adnan matrix estimated']);
    stain_matrix = EstStainUsingMacenko( im );
    [ ~, H2, E2, ~, ~ ] = Deconvolve( im, stain_matrix', verbose );

    disp(['Color Deconvolution using ref image matrix']);
    stain_matrix = EstStainUsingMacenko( Target );
    [ ~, H3, E3, ~, ~ ] = Deconvolve( im, stain_matrix', verbose );

    [NormMM ] = Macenko(im, Target, 255, 0.15, 1, verbose);
    
    imwrite(NormMM, strcat(path, 'Normalized\',imageName,imageExt));
    imwrite(H, strcat(path, 'HE\',imageName,'_H',imageExt));
    imwrite(E, strcat(path, 'HE\',imageName,'_E',imageExt));
    imwrite(H2, strcat(path, 'HE2\',imageName,'_H',imageExt));
    imwrite(E2, strcat(path, 'HE2\',imageName,'_E',imageExt));
    imwrite(H3, strcat(path, 'HE3\',imageName,'_H',imageExt));
    imwrite(E3, strcat(path, 'HE3\',imageName,'_E',imageExt));
    
%     figure,
%     subplot(231); imshow(H);   title('H (C code)');
%     subplot(232); imshow(E);   title('E (C code)');
%     subplot(233); imshow(im);  title('Bg (C code)');
%     subplot(234); imshow(H2);  title('H (Our Implementation)');
%     subplot(235); imshow(E2);  title('E (Our Implementation)');
%     subplot(236); imshow(NormMM); title('Bg (Our Implementation)');
%     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
end

%%

