SNR = 0:2:60;
close all;

%filtre Dirichlet 

[FoundDiracsLocation1,F] = FRI_proj([100 200],0,2,40,10^(-3),SNR,2,1);
[FoundDiracsLocation2,F] = FRI_proj([100 200],0,2,40,10^(-3),SNR,5,1);
[FoundDiracsLocation3,F] = FRI_proj([100 200],0,2,40,10^(-3),SNR,30,1);
[FoundDiracsLocation4,F] = FRI_proj([100 200],0,2,40,10^(-3),SNR,80,1);

figure
plot(-length(F)/2+1:length(F)/2,F)

figure 
plot(SNR,FoundDiracsLocation1(:,1),SNR,FoundDiracsLocation1(:,2),'-s',...
     SNR,FoundDiracsLocation2(:,1),'--',SNR,FoundDiracsLocation2(:,2),'-d',...
     SNR,FoundDiracsLocation3(:,1),'-.',SNR,FoundDiracsLocation3(:,2),'->',...
     SNR,FoundDiracsLocation4(:,1),':',SNR,FoundDiracsLocation4(:,2),'-<');

 %filtre Gaussien

[FoundDiracsLocation1,F] = FRI_proj([100 200],0,1,40,10^(-3),SNR,2,1);
[FoundDiracsLocation2,F] = FRI_proj([100 200],0,1,40,10^(-3),SNR,5,1);
[FoundDiracsLocation3,F] = FRI_proj([100 200],0,1,40,10^(-3),SNR,30,1);
[FoundDiracsLocation4,F] = FRI_proj([100 200],0,1,40,10^(-3),SNR,80,1);


figure 
plot(SNR,FoundDiracsLocation1(:,1),SNR,FoundDiracsLocation1(:,2),'-s',...
     SNR,FoundDiracsLocation2(:,1),'--',SNR,FoundDiracsLocation2(:,2),'-d',...
     SNR,FoundDiracsLocation3(:,1),'-.',SNR,FoundDiracsLocation3(:,2),'->',...
     SNR,FoundDiracsLocation4(:,1),':',SNR,FoundDiracsLocation4(:,2),'-<');

