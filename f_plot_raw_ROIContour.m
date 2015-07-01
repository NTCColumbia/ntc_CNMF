function [weigh_image]=f_plot_raw_ROIContour( A,d1,d2,numberLabel,contourColor,ind_neur )
%F_PLOTROICONTOUR Summary of this function goes here
%   Detailed explanation goes here

nr = size(A,2);   % ROI number
centroid = zeros(nr,2);
if nargin < 6
    ind_neur = ones(nr,1);
    layerNumber = 1;
else
    switchLayer = find(ind_neur==2,1);
    layerNumber = 2;
end

%% plot contour itself

for idx=1:nr
    B=full(reshape(A(:,idx),d1,d2)); 
    [cen_y cen_x]=find_centroid(B);
    centroid(idx,:)=[cen_x cen_y];
    %plot(B{1}(:,2), B{1}(:,1),'Linewidth',2, 'color', contourColor(ind_neur(idx),:));
    %p=patch(B{1}(:,2), B{1}(:,1), contourColor(ind_neur(idx),:));
    %set(p,'FaceAlpha', 0.3);
    %set(p,'edgeColor', [1 0 0]); 
    %set(gcf, 'Renderer', 'OpenGL');
end





%% plot contour with weighted color map
figure;
temp=sum(A,2);
temp=reshape(temp, d1, d2);
weigh_image=temp;
imagesc(temp);
if numberLabel==1 & layerNumber == 1
    for idx=1:nr
        text(centroid(idx,1), centroid(idx,2), num2str(idx),'fontsize',10);
    end
end

if numberLabel==1 & layerNumber == 2
    for idx=1:switchLayer-1
        text(centroid(idx,1), centroid(idx,2), num2str(idx),'fontsize',10);
    end
    for idx=switchLayer:nr
        text(centroid(idx,1), centroid(idx,2), ['b' num2str(idx-switchLayer+1)],'fontsize',10);
%        text(centroid(idx,1), centroid(idx,2), [num2str(idx)],'fontsize',10);
    end
end

set(gca,'YDir','reverse'); 
axis([0 d2 0 d1]);
rectangle('Position', [0, 0, d2, d1], 'linewidth', 1);
box;
% hold on;
% for idx=1:nr
%     B=bwboundaries(full(reshape(A(:,idx),d1,d2))); 
%     centroid(idx,:)=[(max(B{1}(:,2))+min(B{1}(:,2)))/2 (max(B{1}(:,1))+min(B{1}(:,1)))/2];
%     if ind_neur(idx)==2
%         plot(B{1}(:,2), B{1}(:,1),'Linewidth',1.5, 'color', [0 0 0]);
%     else
%         plot(B{1}(:,2), B{1}(:,1),'Linewidth',1.5, 'color', [1 1 1]);        
%     end
% end
end

