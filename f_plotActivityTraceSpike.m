function f_plotActivityTraceSpike( signal_raw, signal_filtered, signal_inferred, signal_spike, frameRate, normalization, separation, labelling )
%PLOTACTIVITYTRACE Summary of this function goes here
%   Detailed explanation goes here
t=1/frameRate*(1:size(signal_raw,2));
index=find(signal_spike==0);
signal_spike(index)=NaN;

figure;
hold on;

if normalization ~=0
    separation=1;
    for idx=1:size(signal_raw,1)
        plot(t,signal_raw(idx,:)/max(signal_raw(idx,:))-idx*separation,'color',[0.8 0.8 0.8],'linewidth',2);
%        plot(t,signal_filtered(idx,:)/max(signal_filtered(idx,:))-idx*separation,'color',[1 0.7 0.7],'linewidth',1.5);
        plot(t,signal_inferred(idx,:)/max(signal_inferred(idx,:))-idx*separation,'color',[0 0 1],'linewidth',1);
    end
    ylabel('Normalized \DeltaF/F');
else
    for idx=1:size(signal_raw,1)
        plot(t,signal_raw(idx,:)-idx*separation,'color',[0.8 0.8 0.8],'linewidth',2);
%        plot(t,signal_filtered(idx,:)-idx*separation,'color',[1 0.7 0.7],'linewidth',1.5);
        plot(t,signal_inferred(idx,:)-idx*separation,'color',[0 0 1],'linewidth',1);
    end    
    ylabel('\DeltaF/F');
end

for idx=1:size(signal_raw,1)
    scatter(t,signal_spike(idx,:)-1-idx*separation-0.2, '.', 'linewidth',2,'MarkerEdgeColor',[0 0 0]);
end    

axis([0 max(t) -separation*(size(signal_raw,1)+2) separation]);
xlabel('Time (s)');

if labelling==1
    for idx=1:size(signal_raw,1)
        text(max(t)*1.05,-idx*separation, num2str(idx));
    end
end

end

