close all;

SNR = -10:2:14;
[FoundDiracsLocation_5,FoundDiracsWeights_5,sf]=weights_locs([100 140 300],[1 1 1],-5,SNR,30);
[FoundDiracsLocation_1,FoundDiracsWeights_1,sf]=weights_locs([100 140 300],[1 1 1],-1,SNR,30);
[FoundDiracsLocation0,FoundDiracsWeights0,sf]=weights_locs([100 140 300],[1 1 1],0,SNR,30);
[FoundDiracsLocation1,FoundDiracsWeights1,sf]=weights_locs([100 140 300],[1 1 1],1,SNR,30);
[FoundDiracsLocation5,FoundDiracsWeights5,sf]=weights_locs([100 140 300],[1 1 1],5,SNR,30);


figure 
plot(SNR,FoundDiracsLocation_5(:,1),SNR,FoundDiracsLocation_1(:,1),'--',SNR,FoundDiracsLocation0(:,1),'-d',SNR,FoundDiracsLocation1(:,1),'->',SNR,FoundDiracsLocation5(:,1),'-<');

figure 
plot(SNR,FoundDiracsLocation_5(:,2),SNR,FoundDiracsLocation_1(:,2),'--',SNR,FoundDiracsLocation0(:,2),'-d',SNR,FoundDiracsLocation1(:,2),'->',SNR,FoundDiracsLocation5(:,2),'-<');

figure 
plot(SNR,FoundDiracsLocation_5(:,3),SNR,FoundDiracsLocation_1(:,3),'--',SNR,FoundDiracsLocation0(:,3),'-d',SNR,FoundDiracsLocation1(:,3),'->',SNR,FoundDiracsLocation5(:,3),'-<');

figure 
plot(SNR,FoundDiracsWeights_5(:,1),SNR,FoundDiracsWeights_1(:,1),'--',SNR,FoundDiracsWeights0(:,1),'-d',SNR,FoundDiracsWeights1(:,1),'->',SNR,FoundDiracsWeights5(:,1),'-<');

figure 
plot(SNR,FoundDiracsWeights_5(:,2),SNR,FoundDiracsWeights_1(:,2),'--',SNR,FoundDiracsWeights0(:,2),'-d',SNR,FoundDiracsWeights1(:,2),'->',SNR,FoundDiracsWeights5(:,2),'-<');

figure 
plot(SNR,FoundDiracsWeights_5(:,3),SNR,FoundDiracsWeights_1(:,3),'--',SNR,FoundDiracsWeights0(:,3),'-d',SNR,FoundDiracsWeights1(:,3),'->',SNR,FoundDiracsWeights5(:,3),'-<');
