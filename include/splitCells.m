%% Source code for:
%% Article title: An automatic method for robust and fast cell detection in bright field images from high-throughput microscopy
%% MS ID        : 7277230959453875 
%% Authors      : Felix Buggenthin, Carsten Marr, Michael Schwarzfischer, Philipp S Hoppe, Oliver Hilsenbeck, Timm Schroeder and Fabian J Theis
%% Journal      : BMC Bioinformatics, September 2013
%% When using this code in your publication, please cite accordingly
%% Copyright (C) 2013 Felix Buggenthin
%   
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

% Parts of this code have been adapted from CellProfiler 1.0
% (https://github.com/CellProfiler/CellProfiler1.0)

function [L,bw2] = splitCells(I_c,bw,MinDiameter,MaxDiameter,maxEcc,docComplement,doMerge)

if docComplement
    I_c = imcomplement(I_c);
    I_c = im2double(I_c);
end

clearborder = 1;

bw = imfill(bw,'holes');

MaximaSuppressionSize = round(MinDiameter/1.5);

MaximaImage = bwulterode(bw);

DistanceTransformedImage = bwdist(~bw);

Overlaid = imimposemin(-DistanceTransformedImage,MaximaImage);

WatershedBoundaries = watershed(Overlaid) > 0;
bw = bw.*WatershedBoundaries;

bw = bwlabel(bw);

bw_woMerge = bwlabel(bw > 0);
bw_woMerge(bw_woMerge > 0) = 1;

if ~doMerge
    bw2 = bw_woMerge;
    L = bwlabel(bw2);
    return
end

%%% Remove bw with no marker in them (this happens occasionally)
%%% This is a very fast way to get pixel indexes for the bw
tmp = regionprops(bw,'PixelIdxList');
for k = 1:length(tmp)
    %%% If there is no maxima in these pixels, exclude object
    if sum(MaximaImage(tmp(k).PixelIdxList)) == 0
        bw(tmp(k).PixelIdxList) = 0;
    end
end


%%% Label the bw
bw = bwlabel(bw);

%%% Merge small bw

if doMerge
    NumberOfbwBeforeMerge = max(bw(:));
    bw = Mergebw(bw,I_c,[MinDiameter MaxDiameter],maxEcc);
    NumberOfbwAfterMerge = max(bw(:));
    NumberOfMergedbw = NumberOfbwBeforeMerge-NumberOfbwAfterMerge;
end

%%% Remove bw along the border of the image (depends on user input)
tmp = bw;
if clearborder
    bw = imclearborder(bw);
end

%%% Relabel the bw
[bw,NumOfbw] = bwlabel(bw > 0);
L = logical(bw);
bw2 = L;
bw2(L ~= 0) = 1;

L = bwlabel(bw2);

end


function bw = Mergebw(bw,I_c,Diameters,MaxEccentricity)

%%% Find the object that we should try to merge with other bw. The object
%%% numbers of these bw are stored in the variable 'MergeIndex'. The bw
%%% that we will try to merge are either the ones that fall below the specified
%%% MinDiameter threshold, or relatively small bw that are above the MaxEccentricity
%%% threshold. These latter bw are likely to be cells where two maxima have been
%%% found and the watershed transform has divided cells into two parts.
MinDiameter = Diameters(1);
MaxDiameter = Diameters(2);

warning('off', 'MATLAB:divideByZero'); %%% Matlab failing atan vs atan2 in regionprops line 672.
props = regionprops(bw,'EquivDiameter','PixelIdxList','Eccentricity');   % Get diameters of the bw
warning('on', 'MATLAB:divideByZero');
EquivDiameters = cat(1,props.EquivDiameter);
Eccentricities = cat(1,props.Eccentricity);
IndexEccentricity = intersect(find(Eccentricities > MaxEccentricity),find(EquivDiameters < (MinDiameter + (MaxDiameter - MinDiameter)/4)));
IndexDiameter = find(EquivDiameters < MinDiameter);
MergeIndex = unique([IndexDiameter;IndexEccentricity]);

% Try to merge until there are no bw left in the 'MergeIndex' list.
[sr,sc] = size(I_c);
counter = round(numel(MergeIndex)/10);
prevcounter = counter;
while ~isempty(MergeIndex)
    counter = round(numel(MergeIndex)/10);
    if counter < prevcounter
        %fprintf('MergeIndex: %i\n',numel(MergeIndex))
        prevcounter = counter;
    end
    
    % Get next object to merge
    CurrentObjectNbr = MergeIndex(1);
    
    %%% Identify neighbors and put their label numbers in a list 'NeighborsNbr'
    %%% Cut a patch so we don't have to work with the entire image
    [r,c] = ind2sub([sr sc],props(CurrentObjectNbr).PixelIdxList);
    rmax = min(sr,max(r) + 3);
    rmin = max(1,min(r) - 3);
    cmax = min(sc,max(c) + 3);
    cmin = max(1,min(c) - 3);
    [xx, yy] = size(bw);
    if(rmax > xx)
        rmax = xx-1;
    end
    if(cmax > yy)
        cmax = yy-1;
    end
    bwPatch = bw(rmin:rmax,cmin:cmax);
    BinaryPatch = double(bwPatch == CurrentObjectNbr);
    GrownBinaryPatch = conv2(BinaryPatch,double(getnhood(strel('disk',2))),'same') > 0;
    Neighbors = bwPatch .*GrownBinaryPatch;
    NeighborsNbr = setdiff(unique(Neighbors(:)),[0 CurrentObjectNbr]);
    
    
    %%% For each neighbor, calculate a set of criteria based on which we decide if to merge.
    %%% Currently, two criteria are used. The first is a Likelihood ratio that indicates whether
    %%% the interface pixels between the object to merge and its neighbor belong to a background
    %%% class or to an object class. The background class and object class are modeled as Gaussian
    %%% distributions with mean and variance estimated from the image. The Likelihood ratio determines
    %%% to which of the distributions the interface voxels most likely belong to. The second criterion
    %%% is the eccentrity of the object resulting from a merge. The more circular, i.e., the lower the
    %%% eccentricity, the better.
    LikelihoodRatio    = zeros(length(NeighborsNbr),1);
    MergedEccentricity = zeros(length(NeighborsNbr),1);
    for j = 1:length(NeighborsNbr)
        
        %%% Get Neigbor number
        CurrentNeighborNbr = NeighborsNbr(j);
        
        %%% Cut patch which contains both original object and the current neighbor
        [r,c] = ind2sub([sr sc],[props(CurrentObjectNbr).PixelIdxList;props(CurrentNeighborNbr).PixelIdxList]);
        rmax = min(sr,max(r) + 3);
        rmin = max(1,min(r) - 3);
        cmax = min(sc,max(c) + 3);
        cmin = max(1,min(c) - 3);
        [xx, yy] = size(bw);
        if(rmax > xx)
            rmax = xx-1;
        end
        if(cmax > yy)
            cmax = yy-1;
        end

        bwPatch = bw(rmin:rmax,cmin:cmax);
        I_cPatch = I_c(rmin:rmax,cmin:cmax);
        
        %%% Identify object interiors, background and interface voxels
        BinaryNeighborPatch      = double(bwPatch == CurrentNeighborNbr);
        BinaryObjectPatch        = double(bwPatch == CurrentObjectNbr);
        GrownBinaryNeighborPatch = conv2(BinaryNeighborPatch,ones(3),'same') > 0;
        GrownBinaryObjectPatch   = conv2(BinaryObjectPatch,ones(3),'same') > 0;
        Interface                = GrownBinaryNeighborPatch.*GrownBinaryObjectPatch;
        Background               = ((GrownBinaryNeighborPatch + GrownBinaryObjectPatch) > 0) - BinaryNeighborPatch - BinaryObjectPatch - Interface;
        WithinObjectIndex        = find(BinaryNeighborPatch + BinaryObjectPatch);
        InterfaceIndex           = find(Interface);
        BackgroundIndex          = find(Background);
        
        %%% Calculate likelihood of the interface belonging to the background or to an object.
        WithinObjectClassMean   = mean(I_cPatch(WithinObjectIndex));
        WithinObjectClassStd    = std(I_cPatch(WithinObjectIndex)) + sqrt(eps);
        BackgroundClassMean     = mean(I_cPatch(BackgroundIndex));
        BackgroundClassStd      = std(I_cPatch(BackgroundIndex)) + sqrt(eps);
        InterfaceMean           = mean(I_cPatch(InterfaceIndex)); %#ok Ignore MLint
        LogLikelihoodObject     = -log(WithinObjectClassStd^2) - (InterfaceMean - WithinObjectClassMean)^2/(2*WithinObjectClassStd^2);
        LogLikelihoodBackground = -log(BackgroundClassStd^2) - (InterfaceMean - BackgroundClassMean)^2/(2*BackgroundClassStd^2);
        LikelihoodRatio(j)      =  LogLikelihoodObject - LogLikelihoodBackground;
        
        %%% Calculate the eccentrity of the object obtained if we merge the current object
        %%% with the current neighbor.
        MergedObject =  double((BinaryNeighborPatch + BinaryObjectPatch + Interface) > 0);
        tmp = regionprops(MergedObject,'Eccentricity');
        MergedEccentricity(j) = tmp(1).Eccentricity;
        
        %%% Get indexes for the interface pixels in original image.
        %%% These indexes are required if we need to merge the object with
        %%% the current neighbor.
        tmp = zeros(size(I_c));
        tmp(rmin:rmax,cmin:cmax) = Interface;
        tmp = regionprops(double(tmp),'PixelIdxList');
        OrigInterfaceIndex{j} = cat(1,tmp.PixelIdxList); %#ok Ignore MLint
    end
    
    %%% Let each feature rank which neighbor to merge with. Then calculate
    %%% a score for each neighbor. If the neighbors is ranked 1st, it will get
    %%% 1 point; 2nd, it will get 2 points; and so on. The lower score the better.
    [ignore,LikelihoodRank]   = sort(LikelihoodRatio,'descend'); %#ok Ignore MLint % The higher the LikelihoodRatio the better
    [ignore,EccentricityRank] = sort(MergedEccentricity,'ascend'); %#ok Ignore MLint % The lower the eccentricity the better
    NeighborScore = zeros(length(NeighborsNbr),1);
    for j = 1:length(NeighborsNbr)
        NeighborScore(j) = find(LikelihoodRank == j) +  find(EccentricityRank == j);
    end
    
    %%% Go through the neighbors, starting with the highest ranked, and merge
    %%% with the first neighbor for which certain basic criteria are fulfilled.
    %%% If no neighbor fulfil the basic criteria, there will be no merge.
    [ignore,TotalRank] = sort(NeighborScore); %#ok Ignore MLint
    for j = 1:length(NeighborsNbr)
        CurrentNeighborNbr = NeighborsNbr(TotalRank(j));
        
        %%% To merge, the interface between bw must be more likely to belong to the object class
        %%% than the background class. The eccentricity of the merged object must also be lower than
        %%% for the original object.
        if LikelihoodRatio(TotalRank(j)) > 0 && MergedEccentricity(TotalRank(j)) < MaxEccentricity%Eccentricities(CurrentObjectNbr)
            
            %%% OK, let's merge!
            %%% Assign the neighbor number to the current object
            bw(props(CurrentObjectNbr).PixelIdxList) = CurrentNeighborNbr;
            
            %%% Assign the neighbor number to the interface pixels between the current object and the neigbor
            bw(OrigInterfaceIndex{TotalRank(j)}) = CurrentNeighborNbr;
            
            %%% Add the pixel indexes to the neigbor index list
            props(CurrentNeighborNbr).PixelIdxList = cat(1,...
                props(CurrentNeighborNbr).PixelIdxList,...
                props(CurrentObjectNbr).PixelIdxList,...
                OrigInterfaceIndex{TotalRank(j)});
            
            %%% Remove the neighbor from the list of bw to be merged (if it's there).
            MergeIndex = setdiff(MergeIndex,CurrentNeighborNbr);
        end
    end
    
    %%% OK, we are done with the current object, let's go to the next
    MergeIndex = MergeIndex(2:end);
end

%%% Finally, relabel the bw
bw = bwlabel(bw > 0);
end

function [xpts,ypts] = getpoints(AxisHandle)

Position = get(AxisHandle,'Position');
FigureHandle = (get(AxisHandle, 'Parent'));
PointHandles = [];
xpts = [];
ypts = [];
NbrOfPoints = 0;
done = 0;
%%% Turns off the CPimagetool function because it interferes with getting
%%% points.
ImageHandle = get(AxisHandle,'children');
set(ImageHandle,'ButtonDownFcn','');

hold on
while ~done;
    
    UserInput = waitforbuttonpress;                            % Wait for user input
    SelectionType = get(FigureHandle,'SelectionType');         % Get information about the last button press
    CharacterType = get(FigureHandle,'CurrentCharacter');      % Get information about the character entered
    
    % Left mouse button was pressed, add a point
    if UserInput == 0 && strcmp(SelectionType,'normal')
        
        % Get the new point and store it
        CurrentPoint  = get(AxisHandle, 'CurrentPoint');
        xpts = [xpts CurrentPoint(2,1)];
        ypts = [ypts CurrentPoint(2,2)];
        NbrOfPoints = NbrOfPoints + 1;
        
        % Plot the new point
        h = plot(CurrentPoint(2,1),CurrentPoint(2,2),'r.');
        set(AxisHandle,'Position',Position)                   % For some reason, Matlab moves the Title text when the first point is plotted, which in turn resizes the image slightly. This line restores the original size of the image
        PointHandles = [PointHandles h];
        
        % If there are any points, and the right mousebutton or the backspace key was pressed, remove a points
    elseif NbrOfPoints > 0 && ((UserInput == 0 && strcmp(SelectionType,'alt')) || (UserInput == 1 && CharacterType == char(8)))   % The ASCII code for backspace is 8
        
        NbrOfPoints = NbrOfPoints - 1;
        xpts = xpts(1:end-1);
        ypts = ypts(1:end-1);
        delete(PointHandles(end));
        PointHandles = PointHandles(1:end-1);
        
   elseif NbrOfPoints >= 3 && UserInput == 1 && CharacterType == char(13)
        done = 1;
        if ~isempty(PointHandles)
            delete(PointHandles)
        end
    end
end
xpts = round(xpts);
ypts = round(ypts);
hold off
set(ImageHandle,'ButtonDownFcn','CPimagetool');
end