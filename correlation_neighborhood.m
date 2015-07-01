function Cn = correlation_neighborhood(Y,d1,d2)

% construct correlation image based on 8-pixel neighborhood
if ndims(Y) == 3
    Y = reshape(Y,d1*d2,size(Y,3));
end

[d,T] = size(Y);
Cn = zeros(d1,d2);
NB = [ 1 0; -1 0; 1 1; 0 1; -1 1; 1 -1; 0 -1; -1 -1];
%figure(2);
for i = 1:d
    loc = [i - d1*(ceil(i/d1)-1), ceil(i/d1)];
    nb_loc = NB + repmat(loc,[8 1]);
    nb_loc(nb_loc(:,1)>d1,:) = [];
    nb_loc(nb_loc(:,1)<1,:) = [];
    nb_loc(nb_loc(:,2)>d2,:) = [];
    nb_loc(nb_loc(:,2)<1,:) = [];
    nb_loc = nb_loc(:,1) + d1*(nb_loc(:,2)-1);
    Cn(loc(1),loc(2)) = mean(corr(Y(i,:)',Y(nb_loc,:)'));
    %if mod(i,d1) == 0
    %    figure(2); imagesc(Cn); title(sprintf('%i',round(i/d1)));
    %end
end
        