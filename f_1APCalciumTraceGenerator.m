function [ temporalWaveform ] = f_1APCalciumTraceGenerator( frameRate, period, gamma )
%F_1APCALCIUMTRACEGENERATOR Summary of this function goes here
%   Detailed explanation goes here
t=0:1/frameRate:period;
temporalWaveform=zeros(1,length(t));

if length(gamma) == 1
    temporalWaveform(1)=1;
    for idx=2:length(t)
        temporalWaveform(idx)=temporalWaveform(idx-1)*gamma;
    end
else if length(gamma) == 2
    temporalWaveform(1)=1;
    temporalWaveform(2)=temporalWaveform(1)*gamma(1);
    for idx=3:length(t)
        temporalWaveform(idx)=temporalWaveform(idx-2)*gamma(2)+temporalWaveform(idx-1)*gamma(1);
    end
    else
        disp('The program currently does not support third order AP system');
    end
end

end

