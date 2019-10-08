close all;

SNR = -10:2:10;

%computation of the weigths and locations with FRI improved technique
[FoundDiracsLocations,FoundDiracsWeights,sf] = weights_locs([100 140 300],[1 1 -0.5],0,SNR,30);

plot(-length(sf)/2+1:length(sf)/2,sf,'LineWidth',2)
%Diracs locations
figure
plot(SNR,FoundDiracsLocations(:,1),'->',SNR,FoundDiracsLocations(:,2),'-<',...
     SNR,FoundDiracsLocations(:,3),'-o');
  
%Diracs weights
figure
plot(SNR,FoundDiracsWeights(:,1),'->',SNR,FoundDiracsWeights(:,2),'-<',...
     SNR,FoundDiracsWeights(:,3),'-o')
 

