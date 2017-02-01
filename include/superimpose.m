function [out] = superimpose(rgb,bin,C)

    [x1 y1 c1] = size(rgb);
    [x2 y2 c2] = size(bin);

    if x1 ~=x2 || y1~=y2
        error('Size of the images are not identical');
    else
        if(c1 == 1)
            out = uint8(zeros(x1, y1, 3));
            out(:,:,2) = uint8(rgb);
        else
            out = rgb;
        end
        for i=1:x1
            for j=1:y1
                if (bin(i,j)>0)
                    out(i,j,1) = C(1)*255;
                    out(i,j,2) = C(2)*255;
                    out(i,j,3) = C(3)*255;
                end
            end
        end
    end
end