function Ath = threshold_components(A,options)

% post processing of spatial components 
% for each component perform the following:  
%   (i)     perform median filtering   
%   (ii)    keep only pixels that contibute up to a level of total energy  
%   (iii)   perform morphological closing  

    defoptions.thr = 0.9999;
    defoptions.se = strel('disk',1);
    defoptions.medw = [3,3];

    if nargin == 1; options = defoptions; end
    if ~isfield(options,'thr'); options.thr = defoptions.thr; end
    if ~isfield(options,'clos_op'); options.se = defoptions.se; end
    if ~isfield(options,'medw'); options.medw = defoptions.medw; end

    [d,nr] = size(A);
    Ath = spalloc(d,nr,nnz(A));
    for i = 1:nr
        A_temp = full(reshape(A(:,i),options.d1,options.d2));
        A_temp = medfilt2(A_temp,options.medw);
        A_temp = A_temp(:);
        if(max(max(A_temp))==0)
            A_temp=full(A(:,i));
            display(['ROI ' num2str(i) ' is small']);
        end
        [temp,ind] = sort(A_temp(:).^2,'ascend'); 
        temp =  cumsum(temp);
        ff = find(temp > (1-options.thr)*temp(end),1,'first');
        BW = zeros(options.d1,options.d2);
        BW(ind(ff:d)) = 1;
        BW = imclose(BW,options.se);
        [L,NUM] = bwlabel(BW,4);
        nrg = zeros(NUM,1);
        for l = 1:NUM
            ff = find(L==l);
            nrg(l) = sum(A(ff,i).^2);
        end
        [~,indm] = max(nrg);
        ff = find(L==indm);
        Ath(ff,i) = A(ff,i);
    end
    
    
end
