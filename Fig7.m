close all;
 
LenF = 2500;

%% filter definition 
 
sigma_w = 40;
h = exp(-pi*(-1249:1250).^2/sigma_w^2);
ffth = abs(fft(h));
M0 = floor(sum((ffth/ffth(1) > 10^(-3)))/2);
F = h';
 
%% first case Dirac pulses located in 100 and 120

DiracsLocations = [100 120];
NbDiracs = length(DiracsLocations);
DiracsWeights = ones(1,NbDiracs);
s = zeros(LenF,1);
s(DiracsLocations) = DiracsWeights;
 
%% filtered signal definition 
sf = real(ifft(fft(s).*fft(F))); 
shat = 1/LenF*fft(sf);
sH = [conj(shat(M0+1:-1:2));shat(1:M0+1)];

%% definition of phi and fequency truncation

phihat = 1/LenF*fft(F);
phiH = [conj(phihat(M0+1:-1:2));phihat(1:M0+1)];
ytrue = sH./phiH;
 
SNR = -10:2:10;
 
Nbreal = 200;
 
T  = 30;
T1 = 3.5; 
T2 = 0.2;

errorl2_denoise     = zeros(1,length(SNR));
errorl2_denoise_cad = zeros(1,length(SNR));

for k = 1:length(SNR)   
 k
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
  kk = 1;
  X = [];
  while (kk <= LenF)
   while (m_filt(kk) <= T1*gamma_estime),
    kk = kk+1;
    if (kk > LenF)
     break;   
    end
   end
      
   if (kk <= LenF)
    Xcur = [];
    while (abs(m_filt(kk)) > T2*gamma_estime),
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
     
     while (m_filt(k1) <= T1*gamma_estime),
      X1 = [X1 k1];
      k1 = k1-1;
     end
    
     X = setdiff(X,X1);
     Xcur = setdiff(Xcur,X1);
     YY = Xcur(Xcur > T1*gamma_estime);
     if length(YY) < 10
      X = setdiff(X,Xcur);
     end
    else
     k1 = kk-1;
     X1 = [];
     while (m_filt(k1) <= T1*gamma_estime),
      X1 = [X1 k1];
      k1 = k1-1;
     end
     X = setdiff(X,X1);
     Xcur = setdiff(Xcur,X1);
     YY = Xcur(Xcur > T1*gamma_estime);
     if (length(YY) < 10)
      X = setdiff(X,Xcur);
     end   
    end  
   end   
  end
     
  Y = (1:LenF).*((abs(m_filt) < T2*gamma_estime)');
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
  
  %Cadzow denoising 
  %we select the central indices 
  indice = zeros(2*T+1,1);
  indice(T+1:2*T+1) = M0+1:M0+1+T;
  indice(1:T)       = M0+1-T:M0;
  
   K     = length(DiracsLocations);
   vhatm = y(indice);
   
   N = 2*T+1;  %we consider a square matrix
   P = T;	   %the matrices have size N-P x P+1. K<=P<=M required.
   
   Nbiter=50;			%number of iterations. 
   Tdenoised=toeplitz(vhatm(P+1:N),vhatm(P+1:-1:1)); %the noisy matrix is the initial estimate
   
   for iter=1:Nbiter
	[U S V]=svd(Tdenoised);
	Tdenoised=U(:,1:K)*S(1:K,1:K)*(V(:,1:K))';	%SVD truncation -> Tdenoised has rank K
	Tdenoised=Toeplitzation(Tdenoised);			% -> Tdenoised is Toeplitz	
   end
   
   %denoised vector
   y_denoised_cad = zeros(2*T+1,1);
   y_denoised_cad(T+1:2*T+1) = Tdenoised(:,1);
   y_denoised_cad(1:T)       = Tdenoised(1,end:-1:2);
   
   
   if  p == 1, 
    errorl2_denoise_cad(k) = 1/sqrt(2*T+1)*norm(ytrue(M0+1-T:M0+1+T)-y_denoised_cad);
     errorl2_denoise(k)    = 1/sqrt(2*T+1)*norm(ytrue(M0+1-T:M0+1+T)-y_denoise(M0+1-T:M0+1+T));
   else
    errorl2_denoise_cad(k) = errorl2_denoise_cad(k)+...
       1/sqrt(2*T+1)*norm(ytrue(M0+1-T:M0+1+T)-y_denoised_cad);
    errorl2_denoise(k) = errorl2_denoise(k)+...
      1/sqrt(2*T+1)*norm(ytrue(M0+1-T:M0+1+T)-y_denoise(M0+1-T:M0+1+T));
   end
 end
end

figure
 X0 = errorl2_denoise/Nbreal;
 X1 = errorl2_denoise_cad/Nbreal;   
 plot(SNR,X0,SNR,X1,'--');
 hold on;
 
%% signal definition
%% second case the Dirac pulses are located at 100, 150, 300, 400
 LenF = 2500;
 DiracsLocations = [100 150 300 400];
 NbDiracs = length(DiracsLocations);
 DiracsWeights = ones(1,NbDiracs);
 s = zeros(LenF,1);
 s(DiracsLocations) = DiracsWeights;

%% definition of phi and fequency truncation

phihat = 1/LenF*fft(F);
phiH = [conj(phihat(M0+1:-1:2));phihat(1:M0+1)];
ytrue = sH./phiH;
 
SNR = -10:2:10;
 
Nbreal = 200;
 
T  = 30;
T1 = 3.5; 
T2 = 0.2;

errorl2_denoise     = zeros(1,length(SNR));
errorl2_denoise_cad = zeros(1,length(SNR));

for k = 1:length(SNR)   
 k
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
  kk = 1;
  X = [];
  while (kk <= LenF)
   while (m_filt(kk) <= T1*gamma_estime),
    kk = kk+1;
    if (kk > LenF)
     break;   
    end
   end
      
   if (kk <= LenF)
    Xcur = [];
    while (abs(m_filt(kk)) > T2*gamma_estime),
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
     
     while (m_filt(k1) <= T1*gamma_estime),
      X1 = [X1 k1];
      k1 = k1-1;
     end
    
     X = setdiff(X,X1);
     Xcur = setdiff(Xcur,X1);
     YY = Xcur(Xcur > T1*gamma_estime);
     if length(YY) < 10
      X = setdiff(X,Xcur);
     end
    else
     k1 = kk-1;
     X1 = [];
     while (m_filt(k1) <= T1*gamma_estime),
      X1 = [X1 k1];
      k1 = k1-1;
     end
     X = setdiff(X,X1);
     Xcur = setdiff(Xcur,X1);
     YY = Xcur(Xcur > T1*gamma_estime);
     if (length(YY) < 10)
      X = setdiff(X,Xcur);
     end   
    end  
   end   
  end
     
  Y = (1:LenF).*((abs(m_filt) < T2*gamma_estime)');
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
  
  %Cadzow denoising 
  %we select the central indices 
  indice = zeros(2*T+1,1);
  indice(T+1:2*T+1) = M0+1:M0+1+T;
  indice(1:T)       = M0+1-T:M0;
  
   K     = length(DiracsLocations);
   vhatm = y(indice);
   
   N = 2*T+1;  %we consider a square matrix
   P = T;	   %the matrices have size N-P x P+1. K<=P<=M required.
   
   Nbiter=50;			%number of iterations. 
   Tdenoised=toeplitz(vhatm(P+1:N),vhatm(P+1:-1:1)); %the noisy matrix is the initial estimate
   
   for iter=1:Nbiter
	[U S V]=svd(Tdenoised);
	Tdenoised=U(:,1:K)*S(1:K,1:K)*(V(:,1:K))';	%SVD truncation -> Tdenoised has rank K
	Tdenoised=Toeplitzation(Tdenoised);			% -> Tdenoised is Toeplitz	
   end
   
   %denoised vector
   y_denoised_cad = zeros(2*T+1,1);
   y_denoised_cad(T+1:2*T+1) = Tdenoised(:,1);
   y_denoised_cad(1:T)       = Tdenoised(1,end:-1:2);
   
   
   if  p == 1, 
    errorl2_denoise_cad(k) = 1/sqrt(2*T+1)*norm(ytrue(M0+1-T:M0+1+T)-y_denoised_cad);
     errorl2_denoise(k)    = 1/sqrt(2*T+1)*norm(ytrue(M0+1-T:M0+1+T)-y_denoise(M0+1-T:M0+1+T));
   else
    errorl2_denoise_cad(k) = errorl2_denoise_cad(k)+...
       1/sqrt(2*T+1)*norm(ytrue(M0+1-T:M0+1+T)-y_denoised_cad);
    errorl2_denoise(k) = errorl2_denoise(k)+...
      1/sqrt(2*T+1)*norm(ytrue(M0+1-T:M0+1+T)-y_denoise(M0+1-T:M0+1+T));
   end
 end
end

X0 = errorl2_denoise/Nbreal;
X1 = errorl2_denoise_cad/Nbreal;   
plot(SNR,X0,'->',SNR,X1,'-<');
hold off;

