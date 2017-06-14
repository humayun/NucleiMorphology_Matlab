function NucleiSegmentationFeature(ImagePath,ImageExt,Method,SplitNuclei,Resolution,GrayLevels,Channels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Humayun Irshad
% NucleiSeg
% SEGMENTATION OF CELL NUCLEI ALGORITHM
% FEATURE EXTRACTION (MORPHOLOGICAL, GRAPH, INTENSITY, TEXTURES)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%addpath('/home/hi41/WS/MatCode/Lib'); 

	%% Source and Destination Path Setting
	if(nargin < 1)
		ImagePath = 'R:/Beck Lab/Figueroa 72 HE/Presentation/Test/';
	end
	SegPath = strcat(ImagePath,'Segmentation/');
	mkdir(SegPath);
	FeaturesPath = strcat(ImagePath,'Features/');
	mkdir(FeaturesPath);

	%% Reading the specific image
	if(nargin < 2)
		ImageExt = '.png';
	end
	srcFiles = dir(strcat(ImagePath,'*',ImageExt));
	
	%% Nuclei Segmentation Method Selection 
	if(nargin<3)
		Method = 2;
	end
		
	%% Enable or disable Nuclei Separation Method
	if(nargin<4)
		SplitNuclei = 1;
	end

	%% Setting image resolution and calculate min and max nuclei sizes (in pixels)
	if(nargin<5)
		Resolution = 0.245;
	end
	MinPixel = 15 / (Resolution * Resolution);
	MaxPixel = 120 / (Resolution * Resolution);

	%% Specify number of graylevels for feature extraction
	if(nargin<6)
		GrayLevels = 64;
	end
	
	%% Specify Channels for Feature Computation
	if(nargin<7)
		Channels = 4;
	end

	for i = 1:length(srcFiles)
		disp(i);

		[~, ImageName, ImageExt] = fileparts(srcFiles(i).name);
		RGB = imread( strcat(ImagePath, ImageName, ImageExt));
		
		% Nuclei Detection & Segmentation
		[nucleiIndexes, BW] = NucleiDetection_HE(RGB, MinPixel, MaxPixel, Method, SplitNuclei);
		imwrite(BW, strcat(SegPath,ImageName,'_Binary.png'));
		
		
		if(length(nucleiIndexes) > 0)
			overlay = superimpose( RGB, BW, [0 1 0]); 
			imwrite(overlay, strcat(SegPath,ImageName, '_Overlay.png'));
			overlay = ImageOverlayPerimeter( RGB, BW); 
			imwrite(overlay, strcat(SegPath,ImageName, '_Perimeter.png'));

			% Compute Featurtes on whole image
			ComputeMorphologicalFeatures( BW, FeaturesPath, ImageName);
			ComputeIntensityFeatures(RGB, BW, Channels, FeaturesPath, ImageName);
			ComputeCMFeatures(RGB, BW, Channels, GrayLevels, FeaturesPath, ImageName);
			ComputeRLFeatures(RGB, BW, Channels, GrayLevels, FeaturesPath, ImageName);
			ComputeGraphFeatures(BW, FeaturesPath, ImageName, SegPath, RGB);

		end
		close all;
	end

% Target Image Stain Estimation
%Target = imread('H:/CodingStuff/Matlab/stain_normalization_toolbox/Ref.png');

%disp(['Color Deconvolution with Standard Stain Matrix ']);
%[ ~, H1, E1, ~, ~ ] = Deconvolve( RGB, [], verbose );

%disp(['Color Deconvolution using Adnan matrix estimated']);
%stain_matrix = EstStainUsingMacenko( RGB );
%[ ~, H2, E2, ~, ~ ] = Deconvolve( RGB, stain_matrix', verbose );

%disp(['Color Deconvolution using ref image matrix']);
%stain_matrix = EstStainUsingMacenko( Target );
%[ ~, H3, E3, ~, ~ ] = Deconvolve( RGB, stain_matrix', verbose );

%RGB  = Macenko(RGB, Target, 255, 0.15, 1, verbose);
%imwrite(RGB, strcat(SegPath,ImageName, '_StainNormalized.png'));


%centroids = nucleiIndexes(:,2:3);
%centroids  = RemoveNucleiNearBorders( RGB, centroids );
% Centroids Filtering 
%priorityThresholdValue = 0;
%centroids = centroids( centroids(:,1) > priorityThresholdValue ,:);
%fig = DisplayDetectedNuclei( RGB, centroids, 2 );
%saveas(fig, strcat(SegPath,ImageName, '_Centroids', ImageExt));


