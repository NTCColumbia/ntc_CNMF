% this program extract df/f

function [ signal_inferred, signal_filtered, signal_raw, signal_inferred_DC, signal_filtered_DC, signal_raw_DC,Y_fres ] = f_signalExtraction_dfof(Y,A,C,b,f,d1,d2,backgroundSubtractionforRaw,baselineRatio)
%MYDATAPROCESSING_SIGNALEXTRACTION Summary of this function goes here
%   Detailed explanation goes here

T = size(C,2);
if ndims(Y) == 3
    Y = reshape(Y,d1*d2,T);
end
nr = size(A,2);     % number of ROIs
nb = size(f,1);     % number of background components
nA = full(sum(A.^2))';  % energy of each row
Yres = Y - A*C - full(b)*f;
Y_r = spdiags(nA,0,nr,nr)\(A'*Yres) + C; 
Y_fres=Y_r-C;                   %residuals from full ROI cf to Raw movie - used to be able to regenrate original raw traces in saved data structure 
signal_inferred=C;                      % spatial weighting on the ROI pixels, background substraction, and denoising
signal_filtered=Y_r;                    % spatial weighting on the ROI pixels, background substraction


%binarize weights of pixels to fraction "weight_thresh" of max (if max
%weight is 0.6, and thresh is .1, set weights > 0.06-->1, < 0.06--> 0
weight_thresh=.2;
A_raw_mask=A;
for idx = 1:nr
    [tempA,ind] = sort(A(:,idx),'ascend');
    temp = cumsum(tempA);
    index_on = find(temp>=weight_thresh*temp(end));
    index_off=find(temp<weight_thresh*temp(end));
    A_raw_mask(ind(index_on),idx)=1;
    A_raw_mask(ind(index_off),idx)=0;
end
[tempI, tempJ, ~]=find(A_raw_mask);
A_contour=sparse(tempI, tempJ, zeros(1,length(tempI))+1, d1*d2, nr);
numberLabel=1;
contourColor=[1 0 0];
f_plot_raw_ROIContour( full(A_contour),d1,d2,numberLabel,contourColor );
title('binarized rawpixel regions','fontsize',16,'fontweight','bold'); drawnow;

nA_contour = full((sum(A_contour.*A_contour)))';  % energy of each row

signal_raw=spdiags(nA_contour,0,nr,nr)\(A_contour'*(Y-full(b)*f));    % no weighting on the ROI pixels, background substraction

signal_inferred_DC=zeros(nr,1);
signal_filtered_DC=zeros(nr,1);
signal_raw_DC=zeros(nr,1);

for idx=1:nr
% inferred signal
    df=(A(:,idx)'*A(:,idx))*C(idx,:);
    temp=sort(df);
    index=find(temp>0);
    temp=temp(index);
    F_res=mean(temp(1:round(baselineRatio*length(temp))));
    signal_inferred(idx,:)=(df-F_res)./(mean(((A(:,idx)'*b)*f))+F_res);
    signal_inferred_DC(idx)=mean(((A(:,idx)'*b)*f))+F_res;
% filtered signal    
    df=(A(:,idx)'*A(:,idx))*signal_filtered(idx,:);
    temp=sort(df);
    index=find(temp>0);
    temp=temp(index);
    F_res=mean(temp(1:round(baselineRatio*length(temp))));    
    signal_filtered(idx,:)=(df-F_res)./(mean(((A(:,idx)'*b)*f))+F_res);
    signal_filtered_DC(idx)=mean(((A(:,idx)'*b)*f))+F_res;
% raw signal: with background substraction    
    df=(A_contour(:,idx)'*A_contour(:,idx))*signal_raw(idx,:);    
    temp=sort(df);
    index=find(temp>0);
    temp=temp(index);
    F_res=mean(temp(1:round(baselineRatio*length(temp))));
    signal_raw(idx,:)=(df-F_res)./(mean(((A_contour(:,idx)'*b)*f))+F_res);
    signal_raw_DC(idx)=mean(((A_contour(:,idx)'*b)*f))+F_res;
end

% really raw data 
    if backgroundSubtractionforRaw==0
        signal_raw=A_contour'*Y;           % no weighting on the ROI pixels, no background substraction
        [signal_raw, signal_raw_DC]=f_calDfof(signal_raw,baselineRatio);
    end
end

