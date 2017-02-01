function [ overlay ] = ImageOverlayPerimeter( im, bw )
%ImageOverlayBorder Summary of this function goes here
%   Overlay image with perimeter of nuclei

    %% Outline the boundary of nuclei
    bw_p = bwperim(bw);
    
    [x y c] = size(im);
    
    if(c == 3)          % for RGB Image
        for i=1:x
            for j=1:y
                if(bw_p(i,j))
                    im(i,j,1) = 0 ;
                    im(i,j,2) = 255;
                    im(i,j,3) = 0;
                end
            end
        end
        overlay = im;
    else            % for Grayscale Image
        r = zeros(x,y);
        b = zeros(x,y);

        %% Create green color contour of nuclei
        rgb = zeros(x,y,3);
        rgb(:,:,1) = r;
        rgb(:,:,2) = bw_p;
        rgb(:,:,3) = b;
        overlay = imfuse(im, rgb, 'blend');
    end
end

