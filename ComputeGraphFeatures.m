function [Features_VD, Features_DT] = ComputeGraphFeatures( bin, FeaturesPath, ImageName, SegPath, rgb)
    
    s = regionprops(logical(bin),'Centroid');
    centroids = cat(1, s.Centroid);
    [n, ~] = size(centroids);
    
    Features_VD = struct('VD_Length',0,'VD_Area',0,'VD_Perimeter',0);
    Features_DT = struct('DT_Area',0,'DT_Perimeter',0);

    if(n>10)
        centroids = double(centroids);
        [V, C] = voronoin(centroids);
        dt = delaunayTriangulation(centroids);

        Voronoi_Features = zeros(size(C,1),3);
        %1st Col=Length, 2nd Col=Area, 3rd Col=Perimeter
        k = 0;
        for j = 1 : size(C ,1)
            ind = C{j}';
            if ind~=1
                k = k + 1;
                Voronoi_Features(k,1) = size(C{j},2);
                [ geom, ~, ~ ] = polygeom(V(ind,1), V(ind,2));
                Voronoi_Features(k,2) = geom(1);
                Voronoi_Features(k,3) = geom(4);
            end
        end
        Features_VD.VD_Length = Voronoi_Features(1:k,1);
        Features_VD.VD_Area = Voronoi_Features(1:k,2);
        Features_VD.VD_Perimeter = Voronoi_Features(1:k,3);

        Delaunay_Features = zeros(size(dt,1),2);
        %1st Col. Area, 2nd Col. Perimeter
        k = 0;
        for j = 1 : size(dt ,1)
            ind = dt(j,:)';
            if ind~=1
                k = k + 1;
                [ geom, ~, ~ ] = polygeom(V(ind,1), V(ind,2));
                Delaunay_Features(k,1) = geom(1);
                Delaunay_Features(k,2) = geom(4);
            end
        end
        Features_DT.DT_Area = Delaunay_Features(1:k,1);
        Features_DT.DT_Perimeter = Delaunay_Features(1:k,2);

        if(nargin >= 3)
            struct2csv(Features_VD, strcat(FeaturesPath, ImageName, '_VD.csv'));
            struct2csv(Features_DT, strcat(FeaturesPath, ImageName, '_DT.csv'));
        end
        
        if(nargin == 5)
            [x, y, c] = size(rgb);
            if(c == 4)  % convert 4 channels RGB to 3 Channels RGB image
                Red = rgb(:,:,1);
                Green = rgb(:,:,2);
                Blue = rgb(:,:,3);
                rgb = cat(3, Red, Green, Blue);
            end
            if(c == 1)  % Convert grayscale to RGB image
                Red = uint8(zeros(x,y));
                Blue = uint8(zeros(x,y));
                rgb = cat(3, Red, rgb, Blue);
            end
            
            warning('off', 'Images:initSize:adjustingMag');

            fig = figure('Visible','off');
            imshow(rgb);
            hold on
            voronoi(centroids(:,1), centroids(:,2));
            triplot(dt,'-g')
            hold off
            saveas(fig, strcat(FeaturesPath, ImageName, '_G_DT_VD.tiff'), 'tiff');
            close(fig);

            fig = figure('Visible','off');
            imshow(rgb);
            hold on
            voronoi(centroids(:,1), centroids(:,2));
            hold off
            saveas(fig, strcat(FeaturesPath, ImageName, '_G_VD.tiff'), 'tiff');
            close(fig);

            fig = figure('Visible','off');
            imshow(rgb);
            hold on
            triplot(dt,'-g')
            hold off
            saveas(fig, strcat(FeaturesPath, ImageName, '_G_DT.tiff'), 'tiff');
            close(fig);
        end
    end    
    Features_VD = struct2table(Features_VD);
    Features_DT = struct2table(Features_DT);
end
