function [weigh_image]=f_plot_raw_ROIContour( A,d1,d2,numberLabel )

nr = size(A,2);   % ROI number
centroid = zeros(nr,2);

%% plot contour itself
for idx=1:nr
    B=full(reshape(A(:,idx),d1,d2)); 
    [cen_y cen_x]=f_find_centroid(B);
    centroid(idx,:)=[cen_x cen_y];
end

%% plot contour with weighted color map
figure;
temp=sum(A,2);
temp=reshape(temp, d1, d2);
weigh_image=temp;
imagesc(temp);
if numberLabel==1
    for idx=1:nr
        text(centroid(idx,1), centroid(idx,2), num2str(idx),'fontsize',10);
    end
end

set(gca,'YDir','reverse'); 
axis([0 d2 0 d1]);
rectangle('Position', [0, 0, d2, d1], 'linewidth', 1);
box;
end

