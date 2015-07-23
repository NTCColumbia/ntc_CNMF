% The code extracts results from the data structure where the results from running CNMF algorithm are saved.
% This code is written by Weijian Yang and Darcy S. Peterka

% note:
% signal_inferred: df/f data, spatial weighting on the ROI pixels, unmixing, background substraction, and denoising
% signal_filtered: df/f data, spatial weighting on the ROI pixels, unmixing, background substraction
% signal_raw: df/f data, no spatial weighting on the ROI pixels, whether background substraction is done depends on user input parameter "backgroundSubtractionforRaw" 

%% Data processing
% please load the saved .mat file into Matlab workspace, and input the data structure in the following line
f_structname=YourMovieName_20150723_173922;
im_inf=imfinfo(f_structname.datawrite.movieFileName);
d2=im_inf(1).Width;
d1=im_inf(1).Height;
blah=size(im_inf);
T= blah(1);
d1=f_structname.datawrite.d1;
d2=f_structname.datawrite.d2;
T=f_structname.datawrite.T;
   
% plot ROI contour
cellID=[1 10 23 24];
numberLabel=1;
contourColor=[1 0 0];
weigh_image=f_plotROIContour( f_structname.A(:,cellID),d1,d2,numberLabel,contourColor );

% View results ROI by ROI
%f_view_patches_mod(Y,A,C,b,f,d1,d2);

%% Plot the signals
normalization=1; % if 0, no normalization, otherwise normalize each trace relative to its max value
separation=4;
labelling=1;
f_plotActivityTrace( f_structname.signal_raw(cellID,:), f_structname.signal_filtered(cellID,:), f_structname.signal_inferred(cellID,:), f_structname.datawrite.frameRate, normalization, separation, labelling);
file_view_patches_mod(f_structname.signal_raw,f_structname.signal_filtered,f_structname.signal_inferred,f_structname.A,f_structname.C,f_structname.b,f_structname.f,d1,d2, f_structname.datawrite.movieFileName);
