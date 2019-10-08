function [Non_Detect]=test_detect(DiracsLocations,cas,SNR,par1,T)
  
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
 
   
 %% Creation stream of Dirac
 NbDiracs = length(DiracsLocations);
 DiracsWeights = ones(1,NbDiracs);
 s = zeros(LenF,1);
 s(DiracsLocations) = DiracsWeights;

 %% signal definition 
 
 sf = real(ifft(fft(s).*fft(F)));
  %% definition of phi and fequency truncation

 phihat = 1/LenF*fft(F);
 phiH   = [conj(phihat(M0+1:-1:2));phihat(1:M0+1)];
 
 Non_Detect = zeros(1,length(par1));
  
 Nbreal = 300; 
 angle_TZroot = zeros(NbDiracs,1);
 
 for nb =1:Nbreal
   noise = randn(length(sf),1);
   m = sigmerge(sf,noise,SNR);
   mhat = 1/LenF*fft(m);

   mhat_tronc = zeros(length(mhat),1);
   mhat_tronc(1:M0+1) = mhat(1:M0+1);
   mhat_tronc(end-M0+1:end) = mhat(end-M0+1:end);
  
   %frequency truncation of m
   m_filt = real(LenF*ifft(mhat_tronc));
  
   %estimation of the noise level
   gamma_estime = median(abs(m_filt-median(m_filt)))/0.6745;
   for k = 1:length(par1)
    m_filt = real(LenF*ifft(mhat_tronc));
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
      while (m_filt(kk) > 0.2*gamma_estime),
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
       if (length(YY) < par1(k))
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
       if (length(YY) < par1(k))
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
     XX = DiracLocations_estime(1);
     YY = DiracLocations_estime(2);
     if (YY > DiracsLocations(2)+10)||(YY < DiracsLocations(2)-10)||...
        (XX > DiracsLocations(1)+10)||(XX < DiracsLocations(1)-10)
        Non_Detect(k) = Non_Detect(k)+1;
     end 
    else %Number of Diracs equal 3
     XX = DiracLocations_estime(1);
     YY = DiracLocations_estime(2);
     ZZ = DiracLocations_estime(3);
     if (ZZ > DiracsLocations(3)+10)||(ZZ < DiracsLocations(3)-10)||...
        (YY > DiracsLocations(2)+10)||(YY < DiracsLocations(2)-10)||...
        (XX > DiracsLocations(1)+10)||(XX < DiracsLocations(1)-10)
        Non_Detect(k) = Non_Detect(k)+1;
    end 
   end
  end
 end 
 Non_Detect = Non_Detect/Nbreal;
end 