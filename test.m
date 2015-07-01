%% construct signal
% Choose p=2 or p=1
p=2;                 % auto regression order
gamma1=1.5;          % 0<gamma1+gamma2<1
gamma2=-0.55;        % gamma2<0
gamma1^2+4*gamma2;   % this has to be larger than 0 to avoid oscillation

p=1;               % auto regression order 
gamma1=0.95;
gamma2=0;

t=1:300;
c=zeros(1,length(t));      % initialize calcium concentration
c(1)=0;              % at t=1, c=0
c(2)=1;              % at t=2, c=1 (spike)

for idx=3:100        % solve the autoregressive problem by iteration
    c(idx)=c(idx-2)*gamma2+c(idx-1)*gamma1;
end

Y1=c;
Y2=zeros(1,length(t));
Y2(10:end)=Y1(1:end-9)*2;
Y2(130:end)=Y1(1:end-129);
Y3=zeros(1,length(t));
Y3(20:end)=Y1(1:end-19);
Y3(110:end)=Y1(1:end-109);
Y3(200:end)=Y1(1:end-199);
figure; plot(t,Y1,'r',t,Y2,'b',t,Y3,'g');

%% solve gamma1 and gamma2 by auto regression
% Choose which signal to use
Y=[Y1; Y1*2];    % two pixels   please choose Y1, or Y2, or Y3
%Y=[Y2; Y2*2];    % two pixels   please choose Y1, or Y2, or Y3
%Y=[Y3; Y3*2];    % two pixels   please choose Y1, or Y2, or Y3

np=size(Y,1);    % pixel number
lags=5;          % auto correlation lag = 5
lags=lags+p;
include_noise=0;
Ycl = mat2cell(Y,ones(np,1),size(Y,2));
XC = cell(np,1);
for j = 1:np
    XC{j} = xcov(Ycl{j},lags,'biased');
end
XC = cell2mat(XC); 

gv = zeros(np*lags,1);
if ~include_noise
    g =  XC(:,lags:-1:1);
    lags = lags - p;
end
A = zeros(np*lags,p);
for i = 1:np
    if ~include_noise
        A((i-1)*lags + (1:lags),:) = toeplitz(g(i,p:p+lags-1),g(i,p:-1:1));
    else
        A((i-1)*lags + (1:lags),:) = toeplitz(XC(i,lags+(1:lags)),XC(i,lags+(1:p))) - sn(i)^2*eye(lags,p);
        gv((i-1)*lags + (1:lags)) = XC(i,lags+2:end)';
    end
end
if ~include_noise
    gv = g(:,p+1:end)';
end
%ph = pinv(A)*gv(:);
ph = A\gv(:);
