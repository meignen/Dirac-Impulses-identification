close all;

%singular values computed without denoising
[svd_ampl02,NbDiracs_fin02] = FRI_proj_Diracs([100 150 300 400],0,-10:2:10,30,1,50);

%same but with denoising
[svd_ampl12,NbDiracs_fin12] = FRI_proj_Diracs([100 150 300 400],1,-10:2:10,30,1,50);

%second example without denoising
[svd_ampl000,NbDiracs_fin000] = FRI_proj_Diracs([100 120 300],0,-10:2:10,30,1,50);
[svd_ampl002,NbDiracs_fin002] = FRI_proj_Diracs([100 120 300],0,-10:2:10,70,1,50);

%same but with denoising
[svd_ampl100,NbDiracs_fin100] = FRI_proj_Diracs([100 120 300],1,-10:2:10,30,1,50);
[svd_ampl102,NbDiracs_fin102] = FRI_proj_Diracs([100 120 300],1,-10:2:10,70,1,50);

SNR = -10:2:10;


figure 
%T=30, (100,150,300,400) without denoising
plot(SNR,svd_ampl02(1,:),SNR,svd_ampl02(2,:),'--',SNR,svd_ampl02(3,:),'-.',...
     SNR,svd_ampl02(4,:),':',SNR,svd_ampl02(5,:),'-*');
 
figure 
%T=30, (100,120,300) without denoising
plot(SNR,svd_ampl000(1,:),SNR,svd_ampl000(2,:),'--',SNR,svd_ampl000(3,:),'-.',...
     SNR,svd_ampl000(4,:),':');
 
figure 
%T=70, (100,120,300) without denoising
plot(SNR,svd_ampl002(1,:),SNR,svd_ampl002(2,:),'--',SNR,svd_ampl002(3,:),'-.',...
     SNR,svd_ampl002(4,:),':');
  
  
 %same but with the denoising procedure
 
 figure 
 %T=30, (100,150,300,400)
 plot(SNR,svd_ampl12(1,:),SNR,svd_ampl12(2,:),'--',SNR,svd_ampl12(3,:),'-.',...
      SNR,svd_ampl12(4,:),':',SNR,svd_ampl12(5,:),'-*');

 figure
 %T=30, (100,120,300)
 plot(SNR,svd_ampl100(1,:),SNR,svd_ampl100(2,:),'--',SNR,svd_ampl100(3,:),'-.',...
      SNR,svd_ampl100(4,:),':');
  
 figure 
 %T=70,(100,120,300)
 plot(SNR,svd_ampl102(1,:),SNR,svd_ampl102(2,:),'--',SNR,svd_ampl102(3,:),'-.',...
      SNR,svd_ampl102(4,:),':');
 