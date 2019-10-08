close all;
%We estimate the Dirac pulses location when the Dirac pulses are located in
%100 and 120
SNR = -10:2:10;

[FoundDiracsLocation,var_DiracsLocation, FoundDiracsLocation_Pencil,...
          var_DiracsLocation_Pencil,FoundDiracsLocation_Music,...
          var_DiracsLocation_Music,ndetect,ndetect_p,ndetect_m,sf] = methods_loc([100 120],SNR,30);

plot(-length(sf)/2+1:length(sf)/2,sf,'LineWidth',2)

figure 
%the locations of the first Dirac pulse
X = FoundDiracsLocation(1,:);
Y = FoundDiracsLocation_Pencil(1,:);
Z = FoundDiracsLocation_Music(1,:);

subplot(2,1,1), plot(SNR,X,SNR,Y,'-d',SNR,Z,'->','LineWidth',2,'MarkerSize',10);
hold on;

X = FoundDiracsLocation(2,:);
Y = FoundDiracsLocation_Pencil(2,:);
Z = FoundDiracsLocation_Music(2,:);

%the location of the second Dirac pulse
plot(SNR,X,SNR,Y,'-d',SNR,Z,'->','LineWidth',2,'MarkerSize',10);
hold off; 
subplot(2,1,2), plot(SNR,ndetect,SNR,ndetect_p,'-d',SNR,ndetect_m,'->','LineWidth',2,'MarkerSize',10);

X = var_DiracsLocation(1,:,:);
Y = var_DiracsLocation_Pencil(1,:,:);
Z = var_DiracsLocation_Music(1,:,:);

figure
subplot(2,1,1), plot(SNR,X,SNR,Y,'-d',SNR,Z,'->','LineWidth',2,'MarkerSize',10);

X = var_DiracsLocation(2,:,:);
Y = var_DiracsLocation_Pencil(2,:,:);
Z = var_DiracsLocation_Music(2,:,:);

subplot(2,1,2), plot(SNR,X,SNR,Y,'-d',SNR,Z,'->','LineWidth',2,'MarkerSize',10);

%We estimate the Dirac pulses location when the Dirac pulses are located in
%100 and 140
[FoundDiracsLocation,var_DiracsLocation, FoundDiracsLocation_Pencil,...
          var_DiracsLocation_Pencil,FoundDiracsLocation_Music,...
          var_DiracsLocation_Music,ndetect,ndetect_p,ndetect_m,sf] = methods_loc([100 140],SNR,30);

figure       

plot(-length(sf)/2+1:length(sf)/2,sf,'LineWidth',2)

figure 
%the locations of the first Dirac pulse
X = FoundDiracsLocation(1,:);
Y = FoundDiracsLocation_Pencil(1,:);
Z = FoundDiracsLocation_Music(1,:);

subplot(2,1,1), plot(SNR,X,SNR,Y,'-d',SNR,Z,'->','LineWidth',2,'MarkerSize',10);
hold on;

X = FoundDiracsLocation(2,:);
Y = FoundDiracsLocation_Pencil(2,:);
Z = FoundDiracsLocation_Music(2,:);

%the location of the second Dirac pulse
plot(SNR,X,SNR,Y,'-d',SNR,Z,'->','LineWidth',2,'MarkerSize',10);
hold off; 
subplot(2,1,2), plot(SNR,ndetect,SNR,ndetect_p,'-d',SNR,ndetect_m,'->','LineWidth',2,'MarkerSize',10);

X = var_DiracsLocation(1,:,:);
Y = var_DiracsLocation_Pencil(1,:,:);
Z = var_DiracsLocation_Music(1,:,:);

figure
subplot(2,1,1), plot(SNR,X,SNR,Y,'-d',SNR,Z,'->','LineWidth',2,'MarkerSize',10);

X = var_DiracsLocation(2,:,:);
Y = var_DiracsLocation_Pencil(2,:,:);
Z = var_DiracsLocation_Music(2,:,:);

subplot(2,1,2), plot(SNR,X,SNR,Y,'-d',SNR,Z,'->','LineWidth',2,'MarkerSize',10);

