close all;
SNR = -10:2:10;

[FoundDiracsWeights,FoundDiracsWeights_locknown,FoundDiracsWeights_estime]=test_weights([100 120],SNR,10,30);


plot(SNR,FoundDiracsWeights_locknown(:,1),'-d',SNR,FoundDiracsWeights(:,1),'->',SNR,FoundDiracsWeights_estime(:,1),'-o');

[FoundDiracsWeights,FoundDiracsWeights_locknown,FoundDiracsWeights_estime]=test_weights([100 140],SNR,10,30);

figure
plot(SNR,FoundDiracsWeights_locknown(:,1),'-d',SNR,FoundDiracsWeights(:,1),'->',SNR,FoundDiracsWeights_estime(:,1),'-o');


