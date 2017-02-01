
function [Norm] = Reinhard(Source, Target, verbose)

% StainNormRienhard: Normalize a Source image (RGB) with respect to a
% Target image (RGB)
% 
% This Routine takes two images as input,in which one is the Source image
% and other is the Target one,and uses Reinhard's method (Refrence below)
% to normalize the stain of Source image
% 
%
% Inputs:
% Source         - RGB image that need to be normalized;
% Target         - RGB REFERENCE image 
% verbose        - (optional)Display Results (including Histograms).
%                   Default value = 0, don't display histograms
%
%
% Output:
% Norm           - Normalized RGB image
%
%
% References:
% E. Reinhard, M. Adhikhmin, B. Gooch, and P. Shirley, “Color transfer
% between images,” IEEE Computer Graphics and Applications, vol. 21(5), pp.
%  34–41, 2001.
%
% Copyright (c) 2013, Adnan Khan
% Department of Computer Science,
% University of Warwick, UK.
% 
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


% RGB to LAB colour space conversion for Source/Ref Image
SourceLab = RGB2Lab(double(Source));    
srcL = SourceLab(:, :, 1);
srca = SourceLab(:, :, 2);
srcb = SourceLab(:, :, 3);

% Mean of Source image in Lab Colourspace
msL = mean(srcL(:));              
msa = mean(srca(:));              
msb = mean(srcb(:));

% Standard deviation of Source image in Lab Colourspace
stdsL = std(srcL(:));         
stdsa = std(srca(:));          
stdsb = std(srcb(:));



    
% RGB to LAB colour space conversion for Target Image
TargetLab = RGB2Lab(double(Target));    
tgtL = TargetLab(:, :, 1);
tgta = TargetLab(:, :, 2);
tgtb = TargetLab(:, :, 3);

% Mean of Target image in Lab Colourspace
mtL = mean(tgtL(:)); 
mta = mean(tgta(:)); 
mtb = mean(tgtb(:));

% Standard deviation of Target image in Lab Colourspace
stdtL = std(tgtL(:));          
stdta = std(tgta(:));          
stdtb = std(tgtb(:));

   
NormLab = zeros(size(Source));

% Normalize each channel based on statistics of source and target images
for i = 1 : size(Source,1) 
    for j = 1 : size(Source,2)
        NormLab(i,j,1) = ((SourceLab(i, j, 1) - msL) * (stdtL/stdsL)) + mtL;
        NormLab(i,j,2) = ((SourceLab(i, j, 2) - msa) * (stdta/stdsa)) + mta;
        NormLab(i,j,3) = ((SourceLab(i, j, 3) - msb) * (stdtb/stdsb)) + mtb;
    end
end

% LAB to RGB conversion
Norm = Lab2RGB(NormLab);   

% Display results if verbose mode is true
if verbose
    figure, 
    subplot(131); imshow(Target); title('Reference Image');
    subplot(132); imshow(Source); title('Source Image');
    subplot(133); imshow(uint8(Norm)); title('Normalized Image');
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
end

end

        

    


