function [weigh_image]=f_plotROIContour( A,d1,d2,numberLabel,contourColor )
%F_PLOTROICONTOUR Summary of this function goes here
%   Detailed explanation goes here

nr = size(A,2);   % ROI number
centroid = zeros(nr,2);

%% plot contour itself
figure; hold on;
for idx=1:nr
    B=bwboundaries(full(reshape(A(:,idx),d1,d2))); 
    centroid(idx,:)=[(max(B{1}(:,2))+min(B{1}(:,2)))/2 (max(B{1}(:,1))+min(B{1}(:,1)))/2];
    plot(B{1}(:,2), B{1}(:,1),'Linewidth',2, 'color', contourColor);
    p=patch(B{1}(:,2), B{1}(:,1), contourColor);
    set(p,'FaceAlpha', 0.3);
    set(p,'edgeColor', [1 0 0]); 
    set(gcf, 'Renderer', 'OpenGL');
end

if numberLabel==1
    for idx=1:nr
        text(centroid(idx,1), centroid(idx,2), num2str(idx),'fontsize',10);
    end
end

set(gca,'YDir','reverse'); 
axis([0 d2 0 d1]);
rectangle('Position', [0, 0, d2, d1], 'linewidth', 1);
box;

%% plot contour with weighted color map
figure;
temp=sum(A,2);
temp=reshape(temp, d1, d2);
weigh_image=temp;
imagesc(temp);
daspect([256 200 1]);
colormap('parula');

end

