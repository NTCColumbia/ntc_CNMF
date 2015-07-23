% The code extracts results from CNMF algorithm.
% This code is written by Weijian Yang and Darcy S. Peterka

% note:
% signal_inferred: df/f data, spatial weighting on the ROI pixels, unmixing, background substraction, and denoising
% signal_filtered: df/f data, spatial weighting on the ROI pixels, unmixing, background substraction
% signal_raw: df/f data, no spatial weighting on the ROI pixels, whether background substraction is done depends on user input parameter "backgroundSubtractionforRaw" 

% Set write_data=1 to save extracted traces (+ all parameters used to run the code) to a .mat file with a unique name (filename + date + CNMF.mat)
write_data_out=1; 
if use_merged==1
    Cdat=Cm;
    Adat=Am;
else
    Cdat=C;
    Adat=A;
end

%% Data processing
% plot ROI contour
numberLabel=1;
contourColor=[1 0 0];
weigh_image=f_plotROIContour( Adat,d1,d2,numberLabel,contourColor );

% Extract calcium traces
backgroundSubtractionforRaw=1;              % user input: background substration for raw data?
baselineRatio=0.35;                         % to obtain df/f, what fraction of total values of trace (sorted ascending) are used to determine DC component of ROI
[ signal_inferred, signal_filtered, signal_raw ,signal_inferred_DC, signal_filtered_DC, signal_raw_DC, Y_fres] = f_signalExtraction_dfof(Y,Adat,Cdat,b,f,d1,d2,backgroundSubtractionforRaw,baselineRatio);

if(write_data_out)
    data_writeNMF(datawrite,weigh_image,Adat, Ain, Cdat, Cin, b, f, signal_raw, signal_filtered, signal_inferred,Y_fres,use_merged);
end

%% Plot the signals
normalization=1; % if 0, no normalization, otherwise normalize each trace relative to its max value
separation=0.5;
labelling=1;
f_plotActivityTrace( signal_raw, signal_filtered, signal_inferred, frameRate, normalization, separation, labelling);
f_view_patches_mod(Yr,Adat,Cdat,b,f,d1,d2);

%% Plot the signals of selected ROIs
cellID=[5, 7, 10, 12];
normalization=0; % if 0, no normalization, otherwise normalize each trace relative to its max value
separation=2;
labelling=1;
f_plotActivityTrace( signal_raw(cellID,:), signal_filtered(cellID,:), signal_inferred(cellID,:), frameRate, normalization, separation, labelling);
f_plotROIContour( Adat(:,cellID),d1,d2,numberLabel,contourColor );
