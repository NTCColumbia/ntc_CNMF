function [A,C,nr,merged_ROIs,P] = merge_ROIs(Y_res,A,b,C,f,P)

% merging of spatially overlapping components that have highly correlated tmeporal activity
% The correlation threshold for merging overlapping components is user specified in P.merge_thr (default value 0.85)
% Inputs:
% Y_res:        residual movie after subtracting all found components
% A:            matrix of spatial components
% b:            spatial background
% C:            matrix of temporal components
% f:            temporal background
% P:            parameter struct

% Outputs:
% A:            matrix of new spatial components
% C:            matrix of new temporal components
% nr:           new number of components
% merged_ROIs:  list of old components that were merged

% Written by:
% Eftychios A. Pnevmatikakis, Simons Foundation, 2015

if ~isfield(P,'merge_thr'); thr = 0.85; else thr = P.merge_thr; end     % merging threshold
if ~isfield(P,'max_merg'); mx = 50; else mx = P.max_merg; end           % maximum merging operations

nr = size(A,2);
[d,T] = size(Y_res);
C_corr = corr(full(C(1:nr,:)'));
FF1 = triu(C_corr)>= thr;                           % find graph of strongly correlated temporal components

A_corr = triu(A(:,1:nr)'*A(:,1:nr));                
A_corr(1:nr+1:nr^2) = 0;
FF2 = A_corr > 0;                                   % find graph of overlapping spatial components

FF3 = and(FF1,FF2);                                 % intersect the two graphs
[l,c] = graph_connected_comp(sparse(FF3+FF3'));     % extract connected components
MC = [];
for i = 1:c
    if length(find(l==i))>1
        MC = [MC,(l==i)'];
    end
end

cor = zeros(size(MC,2),1);
for i = 1:length(cor)
    fm = find(MC(:,i));
    for j1 = 1:length(fm)
        for j2 = j1+1:length(fm)
            cor(i) = cor(i) + C_corr(fm(j1),fm(j2));
        end
    end
end

[~,ind] = sort(cor,'descend');
nm = min(length(ind),mx);   % number of merging operations
merged_ROIs = cell(nm,1);
mc = 30;
A_merged = zeros(d,nm);
C_merged = zeros(nm,T);
if strcmpi(P.method,'constrained_foopsi')
    P_merged.gn = cell(nm,1);
    P_merged.b = cell(nm,1);
    P_merged.c1 = cell(nm,1);
    P_merged.neuron_sn = cell(nm,1);
end
for i = 1:nm
    merged_ROIs{i} = find(MC(:,ind(i)));
    nC = sqrt(sum(C(merged_ROIs{i},:).^2,2));
    A_merged(:,i) = sum(A(:,merged_ROIs{i})*spdiags(nC,0,length(nC),length(nC)),2);    
    Y_res = Y_res + A(:,merged_ROIs{i})*C(merged_ROIs{i},:)+b*f;
    [cc,~,Y_res,Ptemp] = update_temporal_components(real(Y_res),A_merged(:,i),b,median(spdiags(nC,0,length(nC),length(nC))\C(merged_ROIs{i},:)),f,P);
    [aa,bb] = update_spatial_components(Y_res,cc,f,A_merged(:,i),P);
    A_merged(:,i) = aa;
    [cc,~,~,Ptemp] = update_temporal_components(Y_res,A_merged(:,i),bb,cc,f,P);
    if strcmpi(P.method,'constrained_foopsi') || strcmpi(P.method,'MCEM_foopsi')
        P_merged.gn{i} = Ptemp.gn{1};
        P_merged.b{i} = Ptemp.b{1};
        P_merged.c1{i} = Ptemp.c1{1};
        P_merged.neuron_sn{i} = Ptemp.neuron_sn{1};
    end
    
    C_merged(i,:) = cc;
    if i < nm
        Y_res = Y_res - A_merged(:,i)*cc;
    end
end

neur_id = unique(cell2mat(merged_ROIs));

A = [A(:,1:nr),A_merged,A(:,nr+1:end)];
C = [C(1:nr,:);C_merged;C(nr+1:end,:)];
A(:,neur_id) = [];
C(neur_id,:) = [];
if strcmpi(P.method,'constrained_foopsi') || strcmpi(P.method,'MCEM_foopsi')
    P.b(neur_id) = [];
    P.b(nr - length(neur_id) + (1:nm)) = P_merged.b;
    P.gn(neur_id) = [];
    P.gn(nr - length(neur_id) + (1:nm)) = P_merged.gn;
    P.c1(neur_id) = [];
    P.c1(nr - length(neur_id) + (1:nm)) = P_merged.c1;
    P.neuron_sn(neur_id) = [];
    P.neuron_sn(nr - length(neur_id) + (1:nm)) = P_merged.neuron_sn;
end
nr = nr - length(neur_id) + nm;
