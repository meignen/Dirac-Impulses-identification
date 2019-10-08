function [FoundDiracsLocation,FoundDiracsWeights,sf] = weights_locs(DiracsLocations,DiracsWeights,cas,SNR,T)
  
 sigma_w = 40;
 prec = 10^(-3);
 %we create a filter with different skewness
 gaussian = @(x) (1/sqrt((2*pi))*exp(-x.^2/2));
 skewedgaussian = ...
 @(x,alpha) 2*sqrt(2*pi)*gaussian(sqrt(2*pi)/sigma_w*x).*...
                  normcdf(alpha*sqrt(2*pi)/sigma_w*x);

 h = skewedgaussian(-1249:1250,cas);
 ffth = abs(fft(h));
 M0 = floor(sum((ffth/ffth(1) > prec))/2);
 F = h';  
 Fc = toeplitz(F',[F(1) (F(2500:-1:2))']);
 LenF = length(F);
 
 par = 0.2;
 if abs(cas) ==5
  par1 = 5;
 end
 if abs(cas) ==1
  par1 = 15;
 end
 if abs(cas) == 0
  par1 = 10;
 end
   
 %% Creation stream of Dirac
 NbDiracs = length(DiracsLocations);
 s = zeros(LenF,1);
 s(DiracsLocations) = DiracsWeights;

 %% signal definition 
 
 sf = real(ifft(fft(s).*fft(F)));
  %% definition of phi and fequency truncation

 phihat = 1/LenF*fft(F);
 phiH   = [conj(phihat(M0+1:-1:2));phihat(1:M0+1)];
 
 FoundDiracsLocation       = zeros(length(SNR),NbDiracs);
 FoundDiracsWeights        = zeros(length(SNR),NbDiracs);
 
 Nbreal = 200; 
 angle_TZroot = zeros(NbDiracs,1);
 nombre = zeros(1,length(SNR));
 nb_detect = zeros(1,length(SNR));
 
 
 for nb =1:Nbreal
  nb
  for k= 1:length(SNR) 
   noise = randn(length(sf),1);
   m = sigmerge(sf,noise,SNR(k));
   mhat = 1/LenF*fft(m);

   mhat_tronc = zeros(length(mhat),1);
   mhat_tronc(1:M0+1) = mhat(1:M0+1);
   mhat_tronc(end-M0+1:end) = mhat(end-M0+1:end);
  
   %frequency truncation of m
   m_filt = real(LenF*ifft(mhat_tronc));
  
   %estimation of the noise level
   gamma_estime = median(abs(m_filt-median(m_filt)))/0.6745;
  
   kk = 1;
   X = [];
   while (kk <= LenF)
    while (abs(m_filt(kk)) <= 3.5*gamma_estime),
     kk = kk+1;
     if (kk > LenF)
      break;   
     end
    end
      
    if (kk <= LenF)
     Xcur = [];
     while (abs(m_filt(kk)) > par*gamma_estime),
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
     
      while (abs(m_filt(k1)) <= 3.5*gamma_estime),
       X1 = [X1 k1];
       k1 = k1-1;
      end
      X = setdiff(X,X1);
      Xcur = setdiff(Xcur,X1);
      YY = Xcur(Xcur > 3.5*gamma_estime);
   
      if (length(YY) < par1)
       X = setdiff(X,Xcur);
      end
     else
      k1 = kk-1;
      X1 = [];
     
      while (abs(m_filt(k1)) <= 3.5*gamma_estime),
       X1 = [X1 k1];
       k1 = k1-1;
      end
      X = setdiff(X,X1);
      Xcur = setdiff(Xcur,X1);
      YY = Xcur(Xcur > 3.5*gamma_estime);
      if (length(YY) < par1)
       X = setdiff(X,Xcur);
      end   
     end  
    end
   end 
   Y = (1:LenF).*((abs(m_filt) < par*gamma_estime)');
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
      
   mH = [conj(mhat(M0+1:-1:2));mhat(1:M0+1)]; 
   
   y  = mH./phiH;
   
   %indices selection  
   indice = zeros(2*T+1,1);
   indice(T+1:2*T+1) = M0+1:M0+1+T;
   indice(1:T) = M0+1-T:M0;
  
   y_t = y(indice);
   ind = NbDiracs+1;
   A = toeplitz(y_t(ind:ind+2*T-NbDiracs),y_t(ind:-1:ind-NbDiracs)); 
   [U,S,W] = svd(A);
   h = W(:,NbDiracs+1);%one considers the last column of V
   
   %basic computation of the angle, works only with no downsampling
   TZroot = roots(flipud(h));
    
   for p=1:NbDiracs,
    if (angle(TZroot(p)) < 0)
     angle_TZroot(p) = angle(TZroot(p))+2*pi;  
    else
     angle_TZroot(p) = angle(TZroot(p));   
    end
   end
   
   DiracLocations_estime = sort(round(angle_TZroot/(2*pi)*LenF)+1);
   
 
   %assume two Diracs
   if (NbDiracs == 2)
   
    n_detect = 0;
    
    for p = 1:NbDiracs,
     if abs(DiracLocations_estime(p)-DiracsLocations(p)) < 10    
      n_detect = n_detect+1;
     end 
    end
    
    if n_detect == NbDiracs
      a2i = Fc(:,[DiracLocations_estime(1) , DiracLocations_estime(2)])\m;
      FoundDiracsLocation(k,:) =  FoundDiracsLocation(k,:)+(abs([DiracLocations_estime(1) DiracLocations_estime(2)] - DiracsLocations));
      FoundDiracsWeights(k,:) = FoundDiracsWeights(k,:)+ abs(a2i'-DiracsWeights);
      nb_detect(k) = nb_detect(k)+1;
    end
     
    
   else %Number of Diracs equal 3
    n_detect = 0;
    
    for p = 1:NbDiracs,
     if abs(DiracLocations_estime(p)-DiracsLocations(p)) < 10    
      n_detect = n_detect+1;
     end
     
    end
    
    if n_detect == NbDiracs
      a2i = Fc(:,[DiracLocations_estime(1) , DiracLocations_estime(2) , DiracLocations_estime(3)])\m;
      FoundDiracsLocation(k,:) = FoundDiracsLocation(k,:)+(abs([DiracLocations_estime(1) DiracLocations_estime(2) DiracLocations_estime(3)] - DiracsLocations));
      FoundDiracsWeights(k,:) = FoundDiracsWeights(k,:)+ abs(a2i'-DiracsWeights);
      nb_detect(k) = nb_detect(k)+1;
    end
   end
  end
 end
 
 for k = 1:length(SNR)
  FoundDiracsLocation(k,:)       = FoundDiracsLocation(k,:)/nb_detect(k);
  FoundDiracsWeights(k,:)        = FoundDiracsWeights(k,:)/nb_detect(k);
 end 