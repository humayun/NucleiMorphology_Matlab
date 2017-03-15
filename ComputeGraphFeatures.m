%% 
% Written by Humayun Irshad (humayun.irshad@gmail.com)
% Feb 2015
%%
function [GraphFeatures_VD, GraphFeatures_DT] = ComputeGraphFeatures( BW, FeaturesPath, ImageName, SegPath, RGB)
    
    s = regionprops(logical(BW),'Centroid');
    centroids = cat(1, s.Centroid);
    [n, ~] = size(centroids);
    
    GraphFeatures_VD = struct('VD_Length',0,'VD_Area',0,'VD_Perimeter',0);
    GraphFeatures_DT = struct('DT_Area',0,'DT_Perimeter',0);

    if(n>6)
        centroids = double(centroids);
        [V, C] = voronoin(centroids);
        dt = delaunayTriangulation(centroids);

        Voronoi_Features = zeros(size(C,1),3);
        %1st Col=Length, 2nd Col=Area, 3rd Col=Perimeter
        k = 0;
        for j = 1 : size(C ,1)
            ind = C{j}';
            if size(ind,1) > 2
                if sum(isinf(V(ind,1)))>0 || sum(isinf(V(ind,2)))>0 
                    continue;
                end
                k = k + 1;
                Voronoi_Features(k,1) = size(C{j},2);
                [ geom, ~, ~ ] = polygeom(V(ind,1), V(ind,2));
                Voronoi_Features(k,2) = geom(1);
                Voronoi_Features(k,3) = geom(4);
            end
        end
        GraphFeatures_VD.VD_Length = Voronoi_Features(1:k,1);
        GraphFeatures_VD.VD_Area = Voronoi_Features(1:k,2);
        GraphFeatures_VD.VD_Perimeter = Voronoi_Features(1:k,3);

        Delaunay_Features = zeros(size(dt,1),2);
        %1st Col. Area, 2nd Col. Perimeter
        k = 0;
        for j = 1 : size(dt ,1)
            ind = dt(j,:)';
            if size(ind,1) == 3
                %if sum(isinf(V(ind,1)))>0 || sum(isinf(V(ind,2)))>0 
                if sum(isinf(dt.Points(ind,1)))>0 || sum(isinf(dt.Points(ind,2)))>0 
                    continue;
                end
                k = k + 1;
                %[ geom, ~, ~ ] = polygeom(V(ind,1), V(ind,2));
                [ geom, ~, ~ ] = polygeom(dt.Points(ind,1), dt.Points(ind,2));
                Delaunay_Features(k,1) = geom(1);
                Delaunay_Features(k,2) = geom(4);
            end
        end
        GraphFeatures_DT.DT_Area = Delaunay_Features(1:k,1);
        GraphFeatures_DT.DT_Perimeter = Delaunay_Features(1:k,2);

        GraphFeatures_VD = struct2table(GraphFeatures_VD);
        GraphFeatures_DT = struct2table(GraphFeatures_DT);
        
        if(nargin >= 3)
            writetable(GraphFeatures_VD, strcat(FeaturesPath, ImageName, '_VD.csv'));
            writetable(GraphFeatures_DT, strcat(FeaturesPath, ImageName, '_DT.csv'));
        end
        
        if(nargin == 5)
            [x, y, c] = size(RGB);
            if(c == 4)  % convert 4 channels RGB to 3 Channels RGB image
                Red = RGB(:,:,1);
                Green = RGB(:,:,2);
                Blue = RGB(:,:,3);
                RGB = cat(3, Red, Green, Blue);
            end
            if(c == 1)  % Convert grayscale to RGB image
                Red = uint8(zeros(x,y));
                Blue = uint8(zeros(x,y));
                RGB = cat(3, Red, RGB, Blue);
            end
            
            warning('off', 'Images:initSize:adjustingMag');

            fig = figure('Visible','off');
            imshow(RGB);
            hold on
            voronoi(centroids(:,1), centroids(:,2));
            triplot(dt,'-g');
            hold off
            saveas(fig, strcat(SegPath, ImageName, '_G_DT_VD.tiff'), 'tiff');
            close(fig);

            fig = figure('Visible','off');
            imshow(RGB);
            hold on
            voronoi(centroids(:,1), centroids(:,2));
            hold off
            saveas(fig, strcat(SegPath, ImageName, '_G_VD.tiff'), 'tiff');
            close(fig);

            fig = figure('Visible','off');
            imshow(RGB);
            hold on
            triplot(dt,'-g');
            hold off
            saveas(fig, strcat(SegPath, ImageName, '_G_DT.tiff'), 'tiff');
            close(fig);
        end
    end    
end
