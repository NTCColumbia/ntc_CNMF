function [ normalPixel ] = f_unsaturatedPixelFinder( Y, saturationValue )
%F_SATURATIONINDICATOR Summary of this function goes here
%   Detailed explanation goes here

% normalization
    [d1,d2,T] = size(Y);
    d = d1*d2;
    Y = reshape(Y,d,T)/saturationValue;

% find those pixel with its value at 95% of the maximum value
    idx = find(Y<0.9);
    Y(idx) = 0;
    Y(setdiff(1:d1*d2*T,idx)) = 1;
    
    temp = sum(Y,2);
    normalPixel = find(temp<T*0.005);
    
end

