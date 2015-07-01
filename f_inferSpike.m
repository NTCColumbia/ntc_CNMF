function dataSpike = f_inferSpike( data, frameRate, method, stdThreshold0, stdThresholdSlope, lowpassCutoff, gamma, temporalWaveformThreshold )
%INFERspike Summary of this function goes here
%   Detailed explanation goes here
period=20;

% calculate the power of original data
intensity=(sum(data.^1,2));
intensity=intensity./max(intensity);

%% derivative method
if strcmp(method, 'derivative')
% setup low pass filter
    options.nlfilt.value=40;		    % filter length
    options.lpass.value=lowpassCutoff;  % [Hz] low pass cut off frequency
    options.timeRes.value=1/frameRate;	% [s] sampling time
% low pass filter
    if lowpassCutoff<0
        dataFiltered=data;
    else
        dataFiltered=lowpass_filter(data, options);
    end
% run derivative        
    dataDif=gradient(dataFiltered,1,2);
    index=find(dataDif<0); 
    dataDif(index)=0;
    dataSpike=zeros(size(data));
% setup new std threshold
    stdThreshold=stdThreshold0-intensity*stdThresholdSlope;
    index=find(stdThreshold<2.0);
    stdThreshold(index)=2.0;
% thresholding
    for idx=1:size(data,1)
        dataDif(idx,:)=dataDif(idx,:)/max(dataDif(idx,:));      % normalization
        index=find(dataDif(idx,:)>(stdThreshold(idx)*std(dataDif(idx,:))+mean((dataDif(idx,:)))));   
        dataSpike(idx,index)=1;
    end
%% foopsi method
else if strcmp(method, 'foopsi')
    dataFoopsi=run_my_fast_oopsi(data, frameRate);
    dataSpike=zeros(size(data));
    for idx=1:size(data,1)
        index=find(dataFoopsi(idx,:)>(stdThreshold0*std(dataFoopsi(idx,:))+mean((dataFoopsi(idx,:)))));   
        dataSpike(idx,index)=1;
    end    
%% temporal matching method
else if strcmp(method, 'temporalMatching')
    newData=zeros(size(data,1), size(data,2)+length(f_1APCalciumTraceGenerator(frameRate, period, gamma{1}))-1);
    newData(:,1:size(data,2))=data;
    dataSpike=zeros(size(data));
    dataMatch=zeros(size(data));
    remainder=zeros(size(newData));
% setup low pass filter
    options.nlfilt.value=40;		    % filter length
    options.lpass.value=lowpassCutoff;  % [Hz] low pass cut off frequency
    options.timeRes.value=1/frameRate;	% [s] sampling time
% low pass filter
    if lowpassCutoff<0
        dataFiltered=data;
    else
        dataFiltered=lowpass_filter(data, options);
    end   
% run derivative        
    dataDif=gradient(dataFiltered,1,2);
    index=find(dataDif<0); 
    dataDif(index)=0;
% setup new std threshold for derivative
    stdThreshold=stdThreshold0-intensity*stdThresholdSlope;
    index=find(stdThreshold<2.0);
    stdThreshold(index)=2.0;

    for idx=1:size(data,1)
% temporal matching    
        [dataMatch(idx,:), remainder(idx,:)]=deconv(newData(idx,:)/max(newData(idx,:)), f_1APCalciumTraceGenerator(frameRate, period, gamma{idx}));
        index=find(dataMatch(idx,:)>temporalWaveformThreshold);
% derivative
        dataDif(idx,:)=dataDif(idx,:)/max(dataDif(idx,:));      % normalization
        index2=find(dataDif(idx,:)>(stdThreshold(idx)*std(dataDif(idx,:))+mean((dataDif(idx,:)))));  

% original data larger than 0 
        index3=find(newData(idx,:)>0);

% get the intersection of the above methods        
        index4=intersect(intersect(index, index2),index3);
        dataSpike(idx,index4)=1;
    end
end    
end

end

