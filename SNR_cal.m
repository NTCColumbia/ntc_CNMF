num=65;
zz=sMC_Single_500um_001_Cropped17_Jun_2015_22_46_20.signal_raw(num,:);
zz=sMC_Single_500um_001_Cropped17_Jun_2015_22_46_20.signal_filtered(num,:);

T=length(zz);

sn = getSn(zz-mean(zz));
(sum((zz-mean(zz)).^2)-sn^2*T)/(sn^2*T)

ha=zeros(size(zz));
for i=2:7999
ha(i)=zz(i)-zz(i-1);
end
figure()
hist(ha,512);
ll=medfilt1(zz,20);
maxdelta=max(zz)-min(ll)
figure()
plot(zz)
[h1 h2]=hist(ha,512);