function [ H, E, Bg ] = Channels2RGB( I, Stain1, Stain2, Stain3, M, verbose )
% Macenkono: Normalize the appearance of an RGB Source Image to the
% Reference RGB image.
% 
%
% Input:
% I         - RGB Source image;
% Stain1        - deconvolved channel 1;
% Stain2        - deconvolved channel 2;
% Stain3        - deconvolved channel ;
% M         - (optAnal) Color Deconvolution Matrix;
%                        (default Ruifrok & Johnston H&E Matrix)
% verbose   - (optional) Display Results or not?
%             Default value = 0, don't display 
%
%
% Output:
% H         - RGB image for First Stain (Usually Hematoxylin Channel)
% E         - RGB image for Second Stain (Usually Eosin Channel)
% Bg        - RGB image for Third Stain (Usually cross product of first two channels.
%
% References:
% None.
%
%
% Adnan Khan
% Department of Computer Science,
% University of Warwick, UK.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%% Sanity check
if  ~exist('Stain1', 'var') || isempty(Stain1)
  errordlg('Must provide 3 stain channels');
  return;
end

if ~exist('Stain2', 'var')  || isempty(Stain2) 
  errordlg('Must provide 3 stain channels');
  return;
end

if ~exist('Stain3', 'var') || isempty(Stain3) 
  errordlg('Must provide 3 stain channels');
  return;
end

[h, w] = size(Stain1);

% All channels must be of same size
if size(Stain1)~=size(Stain2) | size(Stain2)~=size(Stain3), 
    errordlg('All channels must be of same size'); 
    return; 
end


%% Display results or not?
if ~exist('verbose', 'var') || isempty(verbose)
   verbose = 0; 
end

Bg = [];

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

% the intensity of light entering the specimen (see section 2a of [1])
A = 255;



%% Make stain concentration matrix by stacking the 3 channels
% C = double(cat(2,Stain1(:),Stain2(:),Stain3(:)))';
Stain1 = double(Stain1)/255;
Stain2 = double(Stain2)/255;
Stain3 = double(Stain3)/255;
C = cat(2,Stain1(:),Stain2(:),Stain3(:))';

% C = reshape(DCh', 3, []);

if nargout > 1 || verbose
    H = A*exp(-M(:, 1) * C(1, :));
    H = reshape(H', h, w, 3);
    H = uint8(H);
end

if nargout > 2 || verbose
    E = A*exp(-M(:, 2) * C(2, :));
    E = reshape(E', h, w, 3);
    E = uint8(E);
end

if nargout == 3 || verbose
    Bg = A*exp(-M(:, 3) * C(3, :));
    Bg = reshape(Bg', h, w, 3);
    Bg = uint8(Bg);
end

if verbose,
    figure,
    subplot(141); imshow(I);    title('Source');
    subplot(142);  imshow(H);    title('Hamatoxylin'); 
    subplot(143);  imshow(E);    title('Eosin');
    subplot(144);  imshow(Bg);   title('Background');
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
end

end

