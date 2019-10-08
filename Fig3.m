close all;
 
 %% signal definition
 %first case the Dirac pulses are located in  100 and 120
 LenF = 2500;
 DiracsLocations = [100 200];
 NbDiracs = length(DiracsLocations);
 DiracsWeights = ones(1,NbDiracs);
 s = zeros(LenF,1);
 s(DiracsLocations) = DiracsWeights;
 
 T = [20 30 40 50];
 
 %% filter definition 
 
 sigma_w = 40;
 h = exp(-pi*(-1249:1250).^2/sigma_w^2);
 ffth = abs(fft(h));
  
 M0 = floor(sum((ffth/ffth(1) > 10^(-3)))/2);
 F = h';
 
 %% filtered signal definition 
 sf = real(ifft(fft(s).*fft(F))); 
 
 shat = 1/LenF*fft(sf);
 sH = [conj(shat(M0+1:-1:2));shat(1:M0+1)];
 shat_tronc = zeros(length(sf),1);
 shat_tronc(1:M0+1)= shat(1:M0+1);
 shat_tronc(end-M0+1:end)= shat(end-M0+1:end);
 sf_filt = real(LenF*ifft(shat_tronc)); %the reconstructed frequency truncated signal 

 %% definition of phi and fequency truncation

 phihat = 1/LenF*fft(F);
 phiH   = [conj(phihat(M0+1:-1:2));phihat(1:M0+1)];

 ytrue = sH./phiH;
 %add some noise 
 
 SNR = [-10 0];
 
 Nbreal = 50;
 
 T1 = 3:0.5:5;
 T2 = 0.2:0.05:0.6;
 
 for k = 1:length(SNR)
  errorl2_sf_denoise  = zeros(length(T),length(T1),length(T2));
  errorl2_sf          = zeros(1,length(T));
  
  for p=1:Nbreal, 
   noise = randn(length(sf),1);
   m = sigmerge(sf,noise,SNR(k));

   mhat = 1/LenF*fft(m);
 
   mH = [conj(mhat(M0+1:-1:2));mhat(1:M0+1)];
   y  = mH./phiH;
 
   mhat_tronc = zeros(length(mhat),1);
   mhat_tronc(1:M0+1) = mhat(1:M0+1);
   mhat_tronc(end-M0+1:end) = mhat(end-M0+1:end);
  
   m_filt = real(LenF*ifft(mhat_tronc));
  
   %% denoising step
   gamma_estime = median(abs(m_filt-median(m_filt)))/0.6745;
   
   %% filtering oscillations using piecewise cubic polynomial interpolation
   for q=1:length(T1),
    for r=1:length(T2),   
     kk = 1;
     X = [];
     while (kk <= LenF)
      while (m_filt(kk) <= T1(q)*gamma_estime),
       kk = kk+1;
       if (kk > LenF)
        break;   
       end
      end
      
      if (kk <= LenF)
       Xcur = [];
       while (abs(m_filt(kk)) > T2(r)*gamma_estime),
        X = [X kk];
        Xcur = [Xcur kk];
        kk = kk+1;
        if (kk > LenF),
         break;
        end 
       end
  
       %we move backward to remove the points that were wrongly added
       if (kk <= LenF)
        k1 = kk;
        X1 = [];
     
        while (m_filt(k1) <= T1(q)*gamma_estime),
         X1 = [X1 k1];
         k1 = k1-1;
        end
        X = setdiff(X,X1);
        Xcur = setdiff(Xcur,X1);
        YY = Xcur(Xcur > T1(q)*gamma_estime);
        if length(YY) < 10
         X = setdiff(X,Xcur);
        end
       else
        k1 = kk-1;
        X1 = [];
        while (m_filt(k1) <= T1(q)*gamma_estime),
         X1 = [X1 k1];
         k1 = k1-1;
        end
        X = setdiff(X,X1);
        Xcur = setdiff(Xcur,X1);
        YY = Xcur(Xcur > T1(q)*gamma_estime);
        if (length(YY) < 10)
         X = setdiff(X,Xcur);
        end   
       end  
      end   
     end
     
     Y = (1:LenF).*((abs(m_filt) < T2(r)*gamma_estime)');
     Y = Y(Y>0);
     abscissae = [X Y];
     if ismember(LenF,abscissae)
      ordinates =  [m_filt(X)' zeros(1,length(Y))];
     else
      abscissae = [X Y LenF];
      ordinates =  [m_filt(X)' zeros(1,length(Y)) 0];
     end
     m_filt_denoise = pchip(abscissae,ordinates,1:LenF);
     m_filt_denoise = m_filt_denoise(:);
  
     mhat = 1/LenF*fft(m_filt_denoise);
     mH_denoise = [conj(mhat(M0+1:-1:2));mhat(1:M0+1)];
 
     y_denoise  = mH_denoise./phiH;
     if p == 1 
      for qq = 1:length(T),
       errorl2_sf_denoise(qq,q,r) = 1/sqrt(2*T(qq)+1)*norm(ytrue(M0+1-T(qq):M0+1+T(qq))-y_denoise(M0+1-T(qq):M0+1+T(qq)));
      end
     else    
      for qq = 1:length(T),
       errorl2_sf_denoise(qq,q,r) = errorl2_sf_denoise(qq,q,r)+... 
                                  1/sqrt(2*T(qq)+1)*norm(ytrue(M0+1-T(qq):M0+1+T(qq))-y_denoise(M0+1-T(qq):M0+1+T(qq))); 
      end
     end
    end
   end
   if p == 1 
    for qq = 1:length(T),   
     errorl2_sf(qq) = 1/sqrt(2*T(qq)+1)*norm(ytrue(M0+1-T(qq):M0+1+T(qq))-y(M0+1-T(qq):M0+1+T(qq))); 
    end
   else
    for qq = 1:length(T),
     errorl2_sf(qq) = errorl2_sf(qq)+... 
                                  1/sqrt(2*T(qq)+1)*norm(ytrue(M0+1-T(qq):M0+1+T(qq))-y(M0+1-T(qq):M0+1+T(qq))); 
    end    
   end
  end
  figure
  %error before denoising
  errorl2_sf/Nbreal
  for qq =1:length(T) 
   X0 = zeros(length(T1),length(T2));
   X0(:,:)  = errorl2_sf_denoise(qq,:,:)/Nbreal;
   mesh(T2,T1,X0(:,:));
   hold on;
  end
  hold off 
 end
 
 LenF = 2500;
 DiracsLocations = [100 120];
 NbDiracs = length(DiracsLocations);
 DiracsWeights = ones(1,NbDiracs);
 s = zeros(LenF,1);
 s(DiracsLocations) = DiracsWeights;
 
 %% filter definition 
 
 sigma_w = 40;
 h = exp(-pi*(-1249:1250).^2/sigma_w^2);
 ffth = abs(fft(h));
 M0 = floor(sum((ffth/ffth(1) > 10^(-3)))/2);
 F = h';
 
 %% filtered signal definition 
 sf = real(ifft(fft(s).*fft(F))); 
 
 shat = 1/LenF*fft(sf);
 sH = [conj(shat(M0+1:-1:2));shat(1:M0+1)];
 shat_tronc = zeros(length(sf),1);
 shat_tronc(1:M0+1)= shat(1:M0+1);
 shat_tronc(end-M0+1:end)= shat(end-M0+1:end);
 sf_filt = real(LenF*ifft(shat_tronc)); %the reconstructed frequency truncated signal 

 %% definition of phi and fequency truncation

 phihat = 1/LenF*fft(F);
 phiH = [conj(phihat(M0+1:-1:2));phihat(1:M0+1)];
 
 ytrue = sH./phiH;

 %add some noise 
 
 SNR = [-10 0];
 
 Nbreal = 50;
 
 T1 = 3:0.5:5;
 T2 = 0.2:0.05:0.6;
 
 for k = 1:length(SNR)
  errorl2_sf_denoise  = zeros(length(T),length(T1),length(T2));
  errorl2_sf          = zeros(1,length(T));
  for p=1:Nbreal, 
   noise = randn(length(sf),1);
   m = sigmerge(sf,noise,SNR(k));

   mhat = 1/LenF*fft(m);
 
   mH = [conj(mhat(M0+1:-1:2));mhat(1:M0+1)];
   y  = mH./phiH;
 
   mhat_tronc = zeros(length(mhat),1);
   mhat_tronc(1:M0+1) = mhat(1:M0+1);
   mhat_tronc(end-M0+1:end) = mhat(end-M0+1:end);
  
   m_filt = real(LenF*ifft(mhat_tronc));
  
   %% denoising step
   gamma_estime = median(abs(m_filt-median(m_filt)))/0.6745;
   
   %% filtering oscillations using piecewise cubic polynomial interpolation
   for q=1:length(T1),
    for r=1:length(T2),   
     kk = 1;
     X = [];
     while (kk <= LenF)
      while (m_filt(kk) <= T1(q)*gamma_estime),
       kk = kk+1;
       if (kk > LenF)
        break;   
       end
      end
      
      if (kk <= LenF)
       Xcur = [];
       while (abs(m_filt(kk)) > T2(r)*gamma_estime),
        X = [X kk];
        Xcur = [Xcur kk];
        kk = kk+1;
        if (kk > LenF),
         break;
        end 
       end
  
       %we move backward to remove the points that were wrongly added
       if (kk <= LenF)
        k1 = kk;
        X1 = [];
     
        while (m_filt(k1) <= T1(q)*gamma_estime),
         X1 = [X1 k1];
         k1 = k1-1;
        end
        X = setdiff(X,X1);
        Xcur = setdiff(Xcur,X1);
        YY = Xcur(Xcur > T1(q)*gamma_estime);
        if length(YY) < 10
         X = setdiff(X,Xcur);
        end
       else
        k1 = kk-1;
        X1 = [];
        while (m_filt(k1) <= T1(q)*gamma_estime),
         X1 = [X1 k1];
         k1 = k1-1;
        end
        X = setdiff(X,X1);
        Xcur = setdiff(Xcur,X1);
        YY = Xcur(Xcur > T1(q)*gamma_estime);
        if (length(YY) < 10)
         X = setdiff(X,Xcur);
        end   
       end  
      end   
     end
     
     Y = (1:LenF).*((abs(m_filt) < T2(r)*gamma_estime)');
     Y = Y(Y>0);
     abscissae = [X Y];
     if ismember(LenF,abscissae)
      ordinates =  [m_filt(X)' zeros(1,length(Y))];
     else
      abscissae = [X Y LenF];
      ordinates =  [m_filt(X)' zeros(1,length(Y)) 0];
     end
     m_filt_denoise = pchip(abscissae,ordinates,1:LenF);
     m_filt_denoise = m_filt_denoise(:);
  
     mhat = 1/LenF*fft(m_filt_denoise);
     mH_denoise = [conj(mhat(M0+1:-1:2));mhat(1:M0+1)];
 
     y_denoise  = mH_denoise./phiH;
     if p == 1 
      for qq = 1:length(T),
       errorl2_sf_denoise(qq,q,r) = 1/sqrt(2*T(qq)+1)*norm(ytrue(M0+1-T(qq):M0+1+T(qq))-y_denoise(M0+1-T(qq):M0+1+T(qq)));
      end
     else
      for qq = 1:length(T),
       errorl2_sf_denoise(qq,q,r) = errorl2_sf_denoise(qq,q,r)+...
                                1/sqrt(2*T(qq)+1)*norm(ytrue(M0+1-T(qq):M0+1+T(qq))-y_denoise(M0+1-T(qq):M0+1+T(qq)));
      end
     end
    end
   end
   if p == 1 
    for qq = 1:length(T),   
     errorl2_sf(qq) = 1/sqrt(2*T(qq)+1)*norm(ytrue(M0+1-T(qq):M0+1+T(qq))-y(M0+1-T(qq):M0+1+T(qq))); 
    end
   else
    for qq = 1:length(T),
     errorl2_sf(qq) = errorl2_sf(qq)+... 
                                  1/sqrt(2*T(qq)+1)*norm(ytrue(M0+1-T(qq):M0+1+T(qq))-y(M0+1-T(qq):M0+1+T(qq))); 
    end    
   end 
  end
  figure
  errorl2_sf/Nbreal
  for qq = 1:length(T) 
   X0 = zeros(length(T1),length(T2));
   X0(:,:)  = errorl2_sf_denoise(qq,:,:)/Nbreal;
   mesh(T2,T1,X0(:,:));
   hold on;
  end
  hold off;
 end
 
 %% signal definition
 %first case the Dirac pulses are located in  100, 150, 300 and 400
 LenF = 2500;
 DiracsLocations = [100 150 300 400];
 NbDiracs = length(DiracsLocations);
 DiracsWeights = ones(1,NbDiracs);
 s = zeros(LenF,1);
 s(DiracsLocations) = DiracsWeights;
 
 %% filter definition 
 
 sigma_w = 40;
 h = exp(-pi*(-1249:1250).^2/sigma_w^2);
 ffth = abs(fft(h));
  
 M0 = floor(sum((ffth/ffth(1) > 10^(-3)))/2);
 F = h';
 
 %% filtered signal definition 
 sf = real(ifft(fft(s).*fft(F))); 
 
 shat = 1/LenF*fft(sf);
 sH = [conj(shat(M0+1:-1:2));shat(1:M0+1)];
 shat_tronc = zeros(length(sf),1);
 shat_tronc(1:M0+1)= shat(1:M0+1);
 shat_tronc(end-M0+1:end)= shat(end-M0+1:end);
 sf_filt = real(LenF*ifft(shat_tronc)); %the reconstructed frequency truncated signal 

 %% definition of phi and fequency truncation

 phihat = 1/LenF*fft(F);
 phiH = [conj(phihat(M0+1:-1:2));phihat(1:M0+1)];   
 ytrue = sH./phiH;

 %add some noise 
 
 SNR = [-10 0];
 
 Nbreal = 50;
 
 T1 = 3:0.5:5;
 T2 = 0.2:0.05:0.6;
 for k = 1:length(SNR)
  errorl2_sf_denoise  = zeros(length(T),length(T1),length(T2));
  errorl2_sf = zeros(1,length(T));
  
  for p=1:Nbreal, 
   noise = randn(length(sf),1);
   m = sigmerge(sf,noise,SNR(k));

   mhat = 1/LenF*fft(m);
 
   mH = [conj(mhat(M0+1:-1:2));mhat(1:M0+1)];
   y  = mH./phiH;
 
   mhat_tronc = zeros(length(mhat),1);
   mhat_tronc(1:M0+1) = mhat(1:M0+1);
   mhat_tronc(end-M0+1:end) = mhat(end-M0+1:end);
  
   m_filt = real(LenF*ifft(mhat_tronc));
  
   %% denoising step
   gamma_estime = median(abs(m_filt-median(m_filt)))/0.6745;
   
   %% filtering oscillations using piecewise cubic polynomial interpolation
   for q=1:length(T1),
    for r=1:length(T2),   
     kk = 1;
     X = [];
     while (kk <= LenF)
      while (m_filt(kk) <= T1(q)*gamma_estime),
       kk = kk+1;
       if (kk > LenF)
        break;   
       end
      end
      
      if (kk <= LenF)
       Xcur = [];
       while (abs(m_filt(kk)) > T2(r)*gamma_estime),
        X = [X kk];
        Xcur = [Xcur kk];
        kk = kk+1;
        if (kk > LenF),
         break;
        end 
       end
  
       %we move backward to remove the points that were wrongly added
       if (kk <= LenF)
        k1 = kk;
        X1 = [];
     
        while (m_filt(k1) <= T1(q)*gamma_estime),
         X1 = [X1 k1];
         k1 = k1-1;
        end
        X = setdiff(X,X1);
        Xcur = setdiff(Xcur,X1);
        YY = Xcur(Xcur > T1(q)*gamma_estime);
        if length(YY) < 10
         X = setdiff(X,Xcur);
        end
       else
        k1 = kk-1;
        X1 = [];
        while (m_filt(k1) <= T1(q)*gamma_estime),
         X1 = [X1 k1];
         k1 = k1-1;
        end
        X = setdiff(X,X1);
        Xcur = setdiff(Xcur,X1);
        YY = Xcur(Xcur > T1(q)*gamma_estime);
        if (length(YY) < 10)
         X = setdiff(X,Xcur);
        end   
       end  
      end   
     end
     
     Y = (1:LenF).*((abs(m_filt) < T2(r)*gamma_estime)');
     Y = Y(Y>0);
     abscissae = [X Y];
     if ismember(LenF,abscissae)
      ordinates =  [m_filt(X)' zeros(1,length(Y))];
     else
      abscissae = [X Y LenF];
      ordinates =  [m_filt(X)' zeros(1,length(Y)) 0];
     end
     m_filt_denoise = pchip(abscissae,ordinates,1:LenF);
     m_filt_denoise = m_filt_denoise(:);
  
     mhat = 1/LenF*fft(m_filt_denoise);
     mH_denoise = [conj(mhat(M0+1:-1:2));mhat(1:M0+1)];
 
     y_denoise  = mH_denoise./phiH;
     if p == 1 
      for qq = 1:length(T),
       errorl2_sf_denoise(qq,q,r) = 1/sqrt(2*T(qq)+1)*norm(ytrue(M0+1-T(qq):M0+1+T(qq))-y_denoise(M0+1-T(qq):M0+1+T(qq)));
      end
     else
      %errorl2_sf_denoise(q,r) = errorl2_sf_denoise(q,r)+norm(sf_filt-m_filt_denoise); 
      for qq = 1:length(T),
       errorl2_sf_denoise(qq,q,r) = errorl2_sf_denoise(qq,q,r)+...
                                1/sqrt(2*T(qq)+1)*norm(ytrue(M0+1-T(qq):M0+1+T(qq))-y_denoise(M0+1-T(qq):M0+1+T(qq)));
      end
     end  
    end
   end
   if p == 1 
    for qq = 1:length(T),   
     errorl2_sf(qq) = 1/sqrt(2*T(qq)+1)*norm(ytrue(M0+1-T(qq):M0+1+T(qq))-y(M0+1-T(qq):M0+1+T(qq))); 
    end
   else
    for qq = 1:length(T),
     errorl2_sf(qq) = errorl2_sf(qq)+... 
                                  1/sqrt(2*T(qq)+1)*norm(ytrue(M0+1-T(qq):M0+1+T(qq))-y(M0+1-T(qq):M0+1+T(qq))); 
    end    
   end  
  end
  figure
  errorl2_sf/Nbreal
  for qq =1:length(T) 
   X0 = zeros(length(T1),length(T2));
   X0(:,:)  = errorl2_sf_denoise(qq,:,:)/Nbreal;
   mesh(T2,T1,X0(:,:));
   hold on;
  end
  hold off 
 end
  