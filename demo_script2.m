% The code implements a method for simultaneous source extraction and spike inference from large scale calcium imaging movies.
% The algorithm is presented in more detail here: http://arxiv.org/abs/1409.2903
% This code is written by Eftychios A. Pnevmatikakis and Liam Paninski, 
% followed by modification by Weijian Yang and Darcy S. Peterka, 2015

clearvars %-except 
%% user input
%------- Parameters and data save and export -------
write_data_out=1;                   % user input: set to 1 to save out important params as a mat file; set to zero to just display, and not save to external

%------- tiff file information -------
movieFileName = 'Single_150um_024.tif';    % user input: movie file name (can contain the name of folder)
sframe=1;						    % user input: first frame to read
num2read=3000;					    % user input: number of frames to read
int = 1:(num2read-sframe);          % user input: start frame and end frame to be processed (may want to restrict this due to memory issues.  If larger than number of frames, will be set to max frames)
frameRate = 10;                     % user input: movie frame rate [fps]
bitDepth = 16;                      % user input: bit depth of the movie
saturationValue = 2^bitDepth-1;     % user input: saturation pixel value of the movie

%------- ROI information -------
useROIList = 0;                    % user input: user input the centroid of ROI? No: 0, Yes: 1
ROIListFileName = 0;               % user input: if userROIList==1, input xlsx file name of the ROI list (can contain the name of the folder). 
                                   % The file should contain two columns: x, y centroid of the ROI. This ROIListFile is not used if useROIList=0
ROI_nr =100;                       % user input: if userROIList==0, input number of ROIs to be found                                   
ROI_size = 8;                      % user input: maximum size of ROI in pixels (linear one dimension)
ROI_sigma = 3;                     % user input: std of gaussian (size of ROIs,linear one dimension)               

useMultiscale = 0;                 % user input: whether multiscale ROI is used? No: 0, Yes: 1
                                   % If multiscale ROI is used, the following parameters should be input for the smaller scale 
ROI_nr_smallScale = 50;            % user input: number of small ROIs to be found
ROI_size_smallScale = 4;           % user input: maximum size of small ROI in pixels (linear one dimension)
ROI_sigma_smallScale = 1;          % user input: std of gaussian (size of ROIs,linear one dimension)

%------- Calcium dynamics information -------
calciumAutoRegressionOrder = 2;    % user input: (default: 1) auto regression order for calcium fluorescent dynamics

%------- Parameters for spatial component update -------
growthAlgorithm = 'ellipse';       % user input: stepwise dilation growth or single ellipse coverage search
                                   % choice: 'dilation', stepwise dilation growth
                                   % choice: 'ellipse', single ellipse search
maxGrowthIteration = 20;           % user input: (default: 20) if 'dilation' is used, maximum number of iteration
distanceCoverage = 3;              % user input: (default: 3) if 'ellipse' is used, search size for the ROI growth ; for dendrite, use 10
spatialComponentThreshold = 0.8;   % user input: (default: 0.8) threshold to eliminate pixels during spatial component update 
                                   % This value represents the % of cum sum of the weights of the pixels found within the given ROIs search space that are kept. 
                                   % Values closer to 1.0 mean "keep even the weakest pixels".

%------- Method for temporal component update -------
temporalUpdateMethod = 'constrained_foopsi';   % user input: (default: 'constrained_foopsi')
                                               % choice: 'project', uses plain foopsi algorithm
                                               % choice: 'constrained_foopsi', uses constrained foopsi algorithm
                                               % choice: 'noise_constrained', uses lagrangian foopsi algorithm
											   % choice: 'MCEM_foopsi', uses Monte Carlo EM approach and constrained foopsi 
constrained_foopsi_method = 'cvx';             % user input if temporalUpdateMethod == 'constrained_foopsi': (default: 'cvx')
                                               % choice: 'cvx', use the cvx package available from cvxr.com
                                               % choice: 'dual', use dual ascent
                                               % choice: 'lars' uses the least regression algorithm
                                               % choice: 'spgl1' uses the spgl1 package available from math.ucdavis.edu/~mpf/spgl1/  (usually fastest)
gammaSystemBias = 0.96;                           % user input: (default: 1) system bias for gamma factor (i.e. AR coefficient). Current arpfit.m slightly underestimates tau for all neurons.  This term applies a system wide correction. 
                                               % 1 means no bias. gamma calculated from arpfit.m will be multiplied with this gammaSystemBias
                                               
%------- Parameter for ROI merging -------                                               
merge_thr = 0.8;                               % user input: (default: 0.8) threshold for merging ROI.  The higher number, the more alike the temporal components must be to merge ROIs
use_merged=1;                                  % set to 1 to merge related ROIs                                               

% Restore all parameters in a single struct for easy passing to file writing utility.  
% datawrite structure is not used by anything else.
% Need to run "dataProcessing.m" to generate additional signals and save all these parameters  
if (1)
    datawrite.movieFileName = movieFileName;
    datawrite.sframe=sframe;
    datawrite.num2read=num2read;
    datawrite.int = int;
    datawrite.frameRate = frameRate;
    datawrite.saturationValue = saturationValue;
    datawrite.useROIList = useROIList;
    datawrite.ROIListFileName = ROIListFileName;
    datawrite.ROI_nr = ROI_nr;
    datawrite.ROI_size = ROI_size;
    datawrite.ROI_sigma = ROI_sigma;
    datawrite.useMultiscale = useMultiscale;
    datawrite.ROI_nr_smallScale = ROI_nr_smallScale;
    datawrite.ROI_size_smallScale = ROI_size_smallScale;
    datawrite.ROI_sigma_smallScale = ROI_sigma_smallScale;
    datawrite.calciumAutoRegressionOrder = calciumAutoRegressionOrder;
    datawrite.growthAlgorithm=growthAlgorithm;
    datawrite.maxGrowthIteration=maxGrowthIteration;
    datawrite.distanceCoverage=distanceCoverage;
    datawrite.spatialComponentThreshold=spatialComponentThreshold;
    datawrite.temporalUpdateMethod=temporalUpdateMethod;
    datawrite.constrained_foopsi_method=constrained_foopsi_method;
    datawrite.gammaSystemBias=gammaSystemBias;
    datawrite.merge_thr=merge_thr;
    datawrite.use_merged=use_merged;
end

%% load the tiff movie file
Y = double(bigread2(movieFileName,sframe,num2read));    % load image movie using low level file read.  No *compressed* tiffs!!
[d1,d2,T] = size(Y);                   % d1, d2 are dimensions of the movie, and T is the total number of frames 
d = d1*d2;                             % d is the total number of pixels
int=1:T;

datawrite.d1=d1;
datawrite.d2=d2;
datawrite.T=T;

Y_interp = interp_missing_data(Y);      % interpolate missing data (just for pre-processing)
mis_data = find(Y_interp);
Y(mis_data) = Y_interp(mis_data);       % introduce interpolated values for initialization

Y=Y(:,:,int);
T=length(int);
Y_interp=Y_interp(:,int);

%% fast initialization of spatial components using the greedyROI
params.gSiz = ROI_size;           % maximum size of ROI in pixels (linear one dimension)
params.gSig = ROI_sigma;          % std of gaussian (size of ROIs, linear one dimension)

if useROIList == 1                % ROI list is used
    ROIList = xlsread(ROIListFileName);
    nr = size(ROIList,1);         % number of ROI to be found
    [basis, Cin, center, ~] = greedyROI2d_ROIList(Y, nr, params, ROIList);
else if useMultiscale == 0        % ROI list is not used, multiscale ROI is not used
    nr = ROI_nr;                  % number of ROI to be found
    [basis, Cin, center, ~] = greedyROI2d(Y, nr, params);
    else                          % ROI list is not used, but multiscale ROI is used
        nr = ROI_nr;              % number of ROI to be found (large scale)
        params.gSizSmallScale = ROI_size_smallScale;       % maximum size of small ROI in pixels (linear one dimension)
        params.gSigSmallScale = ROI_sigma_smallScale;      % std of gaussian (size of ROIs,linear one dimension)
        nrSmallScale = ROI_nr_smallScale;                  % number of ROI to be found (small scale)
        [basis, Cin, center, ~] = greedyROI2d_multiscale(Y, nr, nrSmallScale, params);
        nr=nr+nrSmallScale;
    end
end

Ain = sparse(reshape(basis,d,nr));  
Cin = Cin';
clear basis;

% display centers of found components
Cn = correlation_neighborhood(Y,d1,d2);
% Cn = mean(Y,3); %correlation_image(Y); %max(Y,[],3); % image statistic
figure;imagesc(Cn);
colormap('bone'); axis equal; axis tight; hold all;
scatter(center(:,2),center(:,1),'mo');
title('Center of ROIs found from initialization algorithm');
    
%% compute estimates of noise for every pixel and a global time constant
ff = find(sum(Ain,2));            % pixels were greedy method found activity 
p = calciumAutoRegressionOrder;   % order of AR system
options.pixels = ff;
Yr = reshape(Y,d,T);
P = arpfit(Yr,p,options);
[bin,fin] = nnmf(max(Yr-Ain*Cin,0),1);
P.interp = Y_interp;
% remove interpolated values
miss_data_int = find(Y_interp);
Yr(miss_data_int) = P.interp(miss_data_int);

%% find saturation pixel
normalPixel = f_unsaturatedPixelFinder( Y, saturationValue );
saturationPixel = setdiff(1:d,normalPixel);

%% update spatial components
P.d1 = d1; P.d2 = d2; P.thr = spatialComponentThreshold;
if strcmpi(growthAlgorithm,'dilation')
    P.maxGrowthIteration = maxGrowthIteration;
    [A,b] = update_spatial_components2(Yr,Cin,fin,Ain,P);
else
    P.dist = distanceCoverage;
    [A,b] = update_spatial_components(Yr,Cin,fin,Ain,P);
end    

ind_empty = find(sum(A)==0);
A(:,ind_empty) = [];
nr = nr - length(ind_empty);

%% update temporal components
P.method = temporalUpdateMethod;
if strcmpi(temporalUpdateMethod,'constrained_foopsi')
    P.p = calciumAutoRegressionOrder;                           % auto regression order
    P.constrained_foopsi_method = constrained_foopsi_method;    % method used for constrained foopsi
    P.saturationPixel = saturationPixel;                        % saturation pixel
    P.gammaSystemBias = gammaSystemBias;                        % system bias for the AR coefficient
    P.restimate_g = 0;                                          % if P.restimate_g=1, the AR coefficient will be re-estimated in constrained_foopsi 
end
if strcmpi(temporalUpdateMethod,'MCEM_foopsi')
    P.p = calciumAutoRegressionOrder;                           % auto regression order
    P.saturationPixel = saturationPixel;                        % saturation pixel
end
[C,f,Y_res,P] = update_temporal_components(Yr,A,b,Cin,fin,P);
datawrite.gn_preMerge=P.gn;

%% merge found components
P.merge_thr = merge_thr;
if use_merged==1
    Am=A;
    Cm=C;
    repeat=use_merged;
    %set Ck to real, in case solver sends back imaginary components
    while repeat
    [Am,Cm,nr,merged_ROIs,Pm] = merge_ROIs(real(Y_res),Am,b,real(Cm),f,P);
    repeat = ~isempty(merged_ROIs);
    disp(nr);
    end
end

%% order ROI
if use_merged==1
    [A_or,C_or,P_or,str] = order_ROIs(Am,Cm,Pm);
else
    [A_or,C_or,P_or,str] = order_ROIs(A,C,P);
end    
[Coor,json_file] = plot_contours(A_or,reshape(P.sn,d1,d2),[],1);

%% view ROI
if use_merged==1
    f_view_patches_mod(Yr,Am,Cm,b,f,d1,d2);
else
    f_view_patches_mod(Yr,A,C,b,f,d1,d2);
end   

%% construct vedio 
param.skip_frame = 1;
param.ind = [1 2 3 4];
param.make_avi = 0;
make_patch_video(Am,Cm,b,f,Yr,d1,d2,param);
