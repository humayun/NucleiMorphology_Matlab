function [ HE ] = AddThirdStainVector( stain_matrix )


% Deconvolve: Deconvolution of an RGB image into its constituent stain
% channels
%
% 
%
% Input:
% stain_matrix [2 x 3] - Stain matrix with R,G,B components for two stains
% 
%
% Output:
% HE           [3 x 3] - Stain matrix with R,G,B components for three
% stains, where the third stain is obtained by the cross product of the
% first tow stains
% 
%
% References:
%  
%  Ruifrok AC, Johnston DA (2001) Quantification of histochemical stainig
%  by color deconvolution. Analytical and Quantitative Cytology and
%  Histology 23: 291–299.
%  Refer also to the G.Landini's C code.
%
% Copyright (c) 2013, Adnan Khan
% Department of Computer Science,
% University of Warwick, UK.
% 
% 
%


if size(stain_matrix,1)<3    
    stain_matrix= [stain_matrix; 0 0 0];
end

cosx = zeros(1,3);
cosy = zeros(1,3);
cosz = zeros(1,3);
len  = zeros(1,3);

for i = 1 : 3
    % normalise vector length     
    len(i)=sqrt(stain_matrix(i,1)*stain_matrix(i,1) + stain_matrix(i,2)*stain_matrix(i,2) + stain_matrix(i,3)*stain_matrix(i,3));
    if (len(i) ~= 0.0)
        cosx(i)= stain_matrix(i,1)/len(i);
        cosy(i)= stain_matrix(i,2)/len(i);
        cosz(i)= stain_matrix(i,3)/len(i);
    end
end
% If stain 3 is not present
if cosx(3)==0.0 && cosy(3)==0.0 && cosz(3)==0.0    
    % approximate its R channel    
    if ((cosx(1)*cosx(1) + cosx(2)*cosx(2))> 1),
        cosx(3)=0.0;
    else
        cosx(3)=sqrt(1.0-(cosx(1)*cosx(1))-(cosx(2)*cosx(2)));
    end
    
    % approximate its G channel
    if ((cosy(1)*cosy(1) + cosy(2)*cosy(2))> 1)
        cosy(3)=0.0;
    else
        cosy(3)=sqrt(1.0-(cosy(1)*cosy(1))-(cosy(2)*cosy(2)));
    end
    
    % approximate its B channel
    if ((cosz(1)*cosz(1) + cosz(2)*cosz(2))> 1)
        cosz(3)=0.0;
    else
        cosz(3)=sqrt(1.0-(cosz(1)*cosz(1))-(cosz(2)*cosz(2)));
    end
end

% normalize vector length for 3rd stain as well
leng= sqrt(cosx(3)*cosx(3) + cosy(3)*cosy(3) + cosz(3)*cosz(3));
cosx(3)= cosx(3)/leng;
cosy(3)= cosy(3)/leng;
cosz(3)= cosz(3)/leng;

HE = [cosx', cosy', cosz']';

end

