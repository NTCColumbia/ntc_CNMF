

% The code extracts results from NMF algorithm.
% This code is written by Weijian Yang and Darcy S. Peterka

% note:
% signal_inferred: df/f data, spatial weighting on the ROI pixels, unmixing, background substraction, and denoising
% signal_filtered: df/f data, spatial weighting on the ROI pixels, unmixing, background substraction
% signal_raw: df/f data, no spatial weighting on the ROI pixels, whether background substraction is done depends on user input parameter "backgroundSubtractionforRaw" 
% signal_spike: detected spike event

%% Data processing

f_structname=aaMC_Single_VS_500um_02229_Jun_2015_15_35_46;
im_inf=imfinfo(f_structname.datawrite.movieFileName);
d2=im_inf(1).Width;
d1=im_inf(1).Height;
blah=size(im_inf);
   T= blah(1);
   
% plot ROI contour
numberLabel=1;
contourColor=[1 0 0];
weigh_image=f_plotROIContour( f_structname.A,d1,d2,numberLabel,contourColor );

% View results ROI by ROI
%f_view_patches_mod(Y,A,C,b,f,d1,d2);
ROIn=size(f_structname.A,2);





%% Plot the signals
normalization=1; % if 0, no normalization, otherwise normalize each trace relative to its max value
separation=2;
labelling=1;

file_view_patches_mod(f_structname.signal_raw,f_structname.signal_filtered,f_structname.signal_inferred,f_structname.A,f_structname.C,f_structname.b,f_structname.f,d1,d2, f_structname.datawrite.movieFileName);
