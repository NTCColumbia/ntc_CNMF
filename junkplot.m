structtmp=sMC_Single_500um_001_Cropped15_Jun_2015_18_12_36;
hold off;
plot(structtmp.signal_filtered(i,:))
hold on
plot(structtmp.signal_inferred(i,:),'r')
plot(structtmp.signal_raw(i,:),'k')
