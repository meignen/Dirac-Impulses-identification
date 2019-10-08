function [FoundDiracsLocation,var_DiracsLocation, FoundDiracsLocation_Pencil,...
          var_DiracsLocation_Pencil,FoundDiracsLocation_Music,...
          var_DiracsLocation_Music,ndetect,ndetect_p,ndetect_m,sf] = methods_loc(DiracsLocations,SNR,T)

 %% Creation stream of Dirac
 LenF = 2500;
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
 
 %% signal definition
 
 sf = real(ifft(fft(s).*fft(F)));
 
 %% definition of phi and fequency truncation

 phihat = 1/LenF*fft(F);
 phiH  = [conj(phihat(M0+1:-1:2));phihat(1:M0+1)]; 
 
 %% Le filtre phihat tronquée
  
 angle_TZroot = zeros(NbDiracs,1);
 FoundDiracsLocation = zeros(NbDiracs,length(SNR));
 var_DiracsLocation  = zeros(NbDiracs,length(SNR));
 
 FoundDiracsLocation_Pencil = zeros(NbDiracs,length(SNR));
 var_DiracsLocation_Pencil  = zeros(NbDiracs,length(SNR));

 FoundDiracsLocation_Music = zeros(NbDiracs,length(SNR));
 var_DiracsLocation_Music  = zeros(NbDiracs,length(SNR));
 
 Nbreal = 400;
 ndetect   = zeros(1,length(SNR));
 ndetect_p = zeros(1,length(SNR));
 ndetect_m = zeros(1,length(SNR));
 
 for nb =1:Nbreal
  for k= 1:length(SNR) 
      
   noise = randn(length(sf),1);
   m = sigmerge(sf,noise,SNR(k));
   mhat = 1/LenF*fft(m);
   mH  = [conj(mhat(M0+1:-1:2));mhat(1:M0+1)];   
   y = mH./phiH;
   
   mhat_tronc = zeros(length(mhat),1);
   mhat_tronc(1:M0+1) = mhat(1:M0+1);
   mhat_tronc(end-M0+1:end) = mhat(end-M0+1:end);
   
   % frequency truncation of m
   m_filt = real(LenF*ifft(mhat_tronc));
    
   % estimation of noise level
   gamma_estime = median(abs(m_filt-median(m_filt)))/0.6745;
  
   %filtering oscillations using piecewise cubic polynomial interpolation
   kk = 1;
   X = [];
   while (kk <= LenF)
    while (m_filt(kk) <= 3.5*gamma_estime),
     kk = kk+1;
     if (kk > LenF)
      break;   
     end
    end
      
    if (kk <= LenF)
     Xcur = [];
     while (abs(m_filt(kk)) > 0.2*gamma_estime),
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
     
      while (m_filt(k1) <= 3.5*gamma_estime),
       X1 = [X1 k1];
       k1 = k1-1;
      end
      X = setdiff(X,X1);
      Xcur = setdiff(Xcur,X1);
      YY = Xcur(Xcur > 3.5*gamma_estime);
      if length(YY) < 10
       X = setdiff(X,Xcur);
      end
     else
      k1 = kk-1;
      X1 = [];
      while (m_filt(k1) <= 3.5*gamma_estime),
       X1 = [X1 k1];
       k1 = k1-1;
      end
      X = setdiff(X,X1);
      Xcur = setdiff(Xcur,X1);
      YY = Xcur(Xcur > 3.5*gamma_estime);
      if (length(YY) < 10)
       X = setdiff(X,Xcur);
      end     
     end  
    end
   end
     
   Y = (1:LenF).*((abs(m_filt) < 0.2*gamma_estime)');
   Y = Y(Y>0);
   abscissae = [X Y];
   if ismember(LenF,abscissae)
    ordinates =  [m_filt(X)' zeros(1,length(Y))];
   else
    abscissae = [X Y LenF];
    ordinates =  [m_filt(X)' zeros(1,length(Y)) 0];
   end
    
   m_filt = pchip(abscissae,ordinates,1:LenF);
   m_filt = m_filt(:);
   mhat = 1/LenF*fft(m_filt);
    
   mH  = [conj(mhat(M0+1:-1:2));mhat(1:M0+1)];   
   y   = mH./phiH;
    
   % indices selection  
   indice = zeros(2*T+1,1);
   indice(T+1:2*T+1) = M0+1:M0+1+T;
   indice(1:T) = M0+1-T:M0;
  
   y_t = y(indice);
         
   % Computation with TLSA
   ind = NbDiracs+1;
   A = toeplitz(y_t(ind:ind+2*T-NbDiracs),y_t(ind:-1:ind-NbDiracs)); 
   [U,S,W] = svd(A);
   h = W(:,NbDiracs+1);%one considers the last column of V
   % basic computation of the angle, works only with no downsampling
   TZroot = roots(flipud(h));
    
   for p=1:NbDiracs,
    if (angle(TZroot(p)) < 0)
     angle_TZroot(p) = angle(TZroot(p))+2*pi;  
    else         
     angle_TZroot(p) = angle(TZroot(p));
    end
   end
   X = sort(round(angle_TZroot/(2*pi)*LenF)+1);
     
   %Computation with matrix pencil
   %M = M0;			%the matrices have size N-P x P+1. K<=P<=M required.
   N = 2*T+1;
   P = T;
   tau = 1;
   K = NbDiracs;
   vhatm = y(indice);
   Y0=zeros(N-P,P);
   Y1=zeros(N-P,P);
	 
   for kk=0:N-P-1
    Y0(kk+1,:)=vhatm((1:P)+kk);
	Y1(kk+1,:)=vhatm((2:P+1)+kk);
   end
     
   [U,S,V] = svd(Y1,0);
   Sp=S(1:K,1:K); Up=U(:,1:K); Vp=V(:,1:K);
   [U,S,V] = svd(Y0,0);
   Y0 = U(:,1:K)*S(1:K,1:K)*V(:,1:K)';
   Zl = Vp*inv(Sp)*Up'*Y0;
   z = eig(Zl);
   [B,I]=sort(abs(z));
   estimtk = round(angle(z(I(end-K+1:end)))/(2*pi)*LenF+1);
   [estimtk,I] = sort(estimtk);

   %%computation with Music
   estimtk1 = round(-rootmusic(vhatm,K)/(2*pi)*LenF+1);
   [estimtk1,I] = sort(estimtk1);
   
   %test whether we have detected a the Dirac pulses
   n_detect = 0;
   n_detect_p = 0; %number of detection with matrix pencil algorithm
   n_detect_m = 0; %number of detection with music algorithm
   
   for p = 1:NbDiracs,
    if abs(X(p)-DiracsLocations(p)) < 10    
     n_detect = n_detect+1;
    end
    
    if abs(estimtk(p)-DiracsLocations(p)) < 10
     n_detect_p = n_detect_p+1;
    end
    
    if abs(estimtk1(p)-DiracsLocations(p)) < 10
     n_detect_m = n_detect_m+1;
    end
    
   end
         
   if (n_detect == NbDiracs) 
    for p=1:NbDiracs
     FoundDiracsLocation(p,k) = FoundDiracsLocation(p,k)+X(p);
     var_DiracsLocation(p,k)  = var_DiracsLocation(p,k)+ abs(X(p)-DiracsLocations(p));
    end
    ndetect(k) = ndetect(k)+1; 
   end
   
   if (n_detect_p == NbDiracs) 
    for p=1:NbDiracs
     FoundDiracsLocation_Pencil(p,k) = FoundDiracsLocation_Pencil(p,k)+estimtk(p);
     var_DiracsLocation_Pencil(p,k)  = var_DiracsLocation_Pencil(p,k)+ abs(estimtk(p)-DiracsLocations(p));
    end
    ndetect_p(k) = ndetect_p(k)+1; 
   end
   
   if (n_detect_m == NbDiracs)
    for p=1:NbDiracs
     FoundDiracsLocation_Music(p,k) = FoundDiracsLocation_Music(p,k)+estimtk1(p);
     var_DiracsLocation_Music(p,k)  = var_DiracsLocation_Music(p,k)+ abs(estimtk1(p)-DiracsLocations(p));
    end
    ndetect_m(k) = ndetect_m(k)+1; 
   end 
   
  end
 end  
 
 for k = 1:length(SNR)   
  FoundDiracsLocation(:,k) = FoundDiracsLocation(:,k)/ndetect(k);
  var_DiracsLocation(:,k)  = var_DiracsLocation(:,k)/ndetect(k);
 end
 
 for k = 1:length(SNR)   
  FoundDiracsLocation_Pencil(:,k) = FoundDiracsLocation_Pencil(:,k)/ndetect_p(k);
  var_DiracsLocation_Pencil(:,k)  = var_DiracsLocation_Pencil(:,k)/ndetect_p(k);
 end
 
 for k = 1:length(SNR)   
  FoundDiracsLocation_Music(:,k) = FoundDiracsLocation_Music(:,k)/ndetect_m(k);
  var_DiracsLocation_Music(:,k)  = var_DiracsLocation_Music(:,k)/ndetect_m(k);
 end
 
 ndetect = ndetect/Nbreal;
 ndetect_p = ndetect_p/Nbreal;
 ndetect_m= ndetect_m/Nbreal;
 
end
