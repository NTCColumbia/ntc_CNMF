function f_plotActivityTrace( signal_raw, signal_filtered, signal_inferred, frameRate, normalization, separation, labelling )
%PLOTACTIVITYTRACE Summary of this function goes here
%   Detailed explanation goes here
t=1/frameRate*(1:size(signal_raw,2));

figure; hold on;

if normalization ~=0
    separation=1.5;
    for idx=1:size(signal_raw,1)
%        plot(t,signal_raw(idx,:)/max(signal_raw(idx,:))-idx*separation,'color',[0.8 0.8 0.8],'linewidth',2);
        plot(t,signal_filtered(idx,:)/max(signal_filtered(idx,:))-idx*separation,'color',[0 0.44 0.74],'linewidth',1);
        plot(t,signal_inferred(idx,:)/max(signal_inferred(idx,:))-idx*separation,'color',[0.84 0.32 0.1],'linewidth',1);
    end
    ylabel('Normalized \DeltaF/F');
else
    for idx=1:size(signal_raw,1)
%        plot(t,signal_raw(idx,:)-idx*separation,'color',[0.8 0.8 0.8],'linewidth',2);
        plot(t,signal_filtered(idx,:)-idx*separation,'color',[1 0.7 0.7],'linewidth',1.5);
        plot(t,signal_inferred(idx,:)-idx*separation,'color',[0 0 1],'linewidth',1);
    end    
    ylabel('\DeltaF/F');
end
 
axis([0 max(t) -separation*(size(signal_raw,1)+2) separation]);
xlabel('Time (s)');

if labelling==1
    for idx=1:size(signal_raw,1)
        text(max(t)*1.05,-idx*separation, num2str(idx));
    end
end

end

