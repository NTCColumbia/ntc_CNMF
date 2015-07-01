function [ dfof, F ] = f_calDfof( originalData,baselineRatio )
%CALDFOF Summary of this function goes here
%   Detailed explanation goes here
    dfof=originalData;
    F=zeros(size(originalData,1),1);
    for idx=1:size(originalData,1)
        temp=sort(originalData(idx,:));
        index=find(temp>0);
        temp=temp(index);
        F(idx)=mean(temp(1:round(baselineRatio*length(temp))));
        dfof(idx,:)=(originalData(idx,:)-F(idx))/F(idx);    
    end    
end

