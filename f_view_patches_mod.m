function f_view_patches_mod(Y,A,C,b,f,d1,d2,ind_neur)

T = size(C,2);
if ndims(Y) == 3
    Y = reshape(Y,d1*d2,T);
end
nr = size(A,2);     % number of ROIs
nb = size(f,1);     % number of background components
nA = full(sum(A.^2))';  % energy of each row
Yres = Y - A*C - full(b)*f;
Y_r = spdiags(nA,0,nr,nr)\(A'*Yres) + C; 
    
%Y_r = spdiags(nA,0,nr,nr)\(A'*(Y-full(b)*f)); 
    
% [~,ind_sor] = sort(max(C,[],2),'descend');
% A = A(:,ind_sor);
% C = C(ind_sor,:);

if nargin < 8
    ind_neur = ones(nr,1);
end

%ind_neur = ind_neur(ind_sor);

figure;
    set(gcf,'Position',[300,300,960,480]);
i=1;
while i <= nr+nb
    subplot(211);
    if i <= nr
        imagesc(reshape(A(:,i),d1,d2)); axis equal; axis tight;
        title({sprintf('Plane %i',ind_neur(i));sprintf('Component %i (press any key to continue)',i)},'fontsize',16,'fontweight','bold'); drawnow; %pause;
    else
        imagesc(reshape(b(:,i-nr),d1,d2)); axis equal; axis tight;
        title('Background component','fontsize',16,'fontweight','bold'); drawnow; 
    end
    subplot(212);
    if i <= nr
%        plot(1:T,Y_r(ind_sor(i),:),1:T,C(i,:)); 
        plot(1:T,Y_r(i,:),1:T,C(i,:)); 
        title(sprintf('Component %i (calcium)',i),'fontsize',16,'fontweight','bold');
        legend('Raw trace (filtered)','Inferred');
        drawnow; 
        a=input('Press ''b'' to go back, or any other key to go forward, followed by ENTER ','s');
        if ~strcmpi(a,'b')
            i = i+1;
        else
            i = max(1,i-1);
        end
    else
        plot(1:T,f(i-nr,:)); title('Background activity','fontsize',16,'fontweight','bold');
        drawnow; 
        a=input('Press ''b'' to go back, or any other key to go forward, followed by ENTER ','s');
        if ~strcmpi(a,'b')
            i = i+1;
        else
            i = max(1,i-1);
        end 
    end
end