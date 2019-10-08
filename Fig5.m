close all;
SNR = -10:2:10;

[svd_ampl10,NbDiracs_fin10] = FRI_proj_Diracs([100 150 300 400],1,SNR,30,1,200);

[svd_ampl100,NbDiracs_fin100] = FRI_proj_Diracs([100 120 300],1,SNR,30,1,200);
[svd_ampl101,NbDiracs_fin101] = FRI_proj_Diracs([100 120 300],1,SNR,70,1,200);


figure 
plot(SNR,NbDiracs_fin10(2,:),SNR,NbDiracs_fin10(3,:),'--');

figure 
plot(SNR,NbDiracs_fin100(2,:),SNR,NbDiracs_fin100(3,:),'--');
hold on
plot(SNR,NbDiracs_fin101(2,:),'-<',SNR,NbDiracs_fin101(3,:),'->');
hold off