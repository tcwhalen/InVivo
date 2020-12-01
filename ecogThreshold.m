function [ecog] = ecogThreshold(ecog,thresh,around)
% Tim Whalen, last edited Nov 2020
% Removes noisy ECoG epochs by thresholding and thresholding

arguments
    ecog (:,1) double % ECoG signal to denoise
    thresh double % magnitude over which ECoG is considered noise, need >0
    around double = 1250 % number of bins around suprathreshold bins which 
        % should also be labeled noise. Default is such that 0.25 seconds
        % removed for 5 kHz signal
end

badBins = abs(ecog)>thresh;
badBins = conv(badBins,ones(2*around,1),'same');
ecog(badBins>=1)=NaN;
end

