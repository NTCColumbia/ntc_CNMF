function [ cen_i cen_j ] = f_find_centroid( image )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
image=image/sum(image(:));
[m,n]=size(image);
[I,J]=ndgrid(1:m,1:n);
centroid=[dot(I(:),image(:)),  dot(J(:),image(:))];  
cen_i=cast(centroid(1),'uint16');
cen_j=cast(centroid(2),'uint16');
end

