function [svd_ampl,NbDiracs_fin] = FRI_proj_Diracs(DiracsLocations,cas,SNR,L,down_samp,Nbreal)
 

%% Creation stream of Dirac
 LenF = 2500;
 NbDiracs = 10;
 DiracsWeights = ones(1,length(DiracsLocations));
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
 phiH   = zeros(length(phihat),1);
 phiH(1:2*M0+1) = [conj(phihat(M0+1:-1:2));phihat(1:M0+1)]; 
 
 %% Dphi
 Dphi = phiH(1:2*M0+1);   
   
 svd_ampl = zeros(NbDiracs+1,length(SNR));  
 NbDiracs_fin = zeros(3,length(SNR));
 
 for nb =1:Nbreal
  for k= 1:length(SNR) 
   noise = randn(length(sf),1);
   m = sigmerge(sf,noise,SNR(k));
   mhat = 1/LenF*fft(m);

   mH = [conj(mhat(M0+1:-1:2));mhat(1:M0+1)];
   
   mhat_tronc = zeros(length(mhat),1);
   mhat_tronc(1:M0+1) = mhat(1:M0+1);
   mhat_tronc(end-M0+1:end) = mhat(end-M0+1:end);
  
   %frequency truncation of m
   m_filt = real(LenF*ifft(mhat_tronc));
   m_filt1 = m_filt;
   
   if (cas == 2)||(cas == 1)
    %estimation of the noise level
    gamma_estime = median(abs(m_filt-median(m_filt)))/0.6745;
  
    %filtering oscillations using piecewise cubic polynomial interpolation
    kk = 1;
    X = [];
    while (kk <= LenF)
     if (kk <= LenF)
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
   end 
   
   y  = mH./Dphi;
   
   indice = zeros(2*L+1,1);
   
   %we select the indices  
   indice(L+1:2*L+1) = M0+1:down_samp:M0+1+L*down_samp;
   indice(1:L) = M0+1-L*down_samp:down_samp:M0+1-down_samp;
  
   y_t = y(indice);
   ind = NbDiracs+1;
   A = toeplitz(y_t(ind:ind+2*L-NbDiracs),y_t(ind:-1:ind-NbDiracs)); 
   %cas
  
   [U,S,W] = svd(A);
   
   %diag(S)
   index2 = (cumsum(diag(S))/sum(diag(S)) > 0.99)'.*(1:NbDiracs+1);
   index1 = (cumsum(diag(S))/sum(diag(S)) > 0.98)'.*(1:NbDiracs+1);
   index0 = (cumsum(diag(S))/sum(diag(S)) > 0.97)'.*(1:NbDiracs+1);
   
   NbDiracs_fin(1,k) = NbDiracs_fin(1,k)+(min(index0(index0> 0))==length(DiracsLocations)); 
   NbDiracs_fin(2,k) = NbDiracs_fin(2,k)+(min(index1(index1 > 0))==length(DiracsLocations)); 
   NbDiracs_fin(3,k) = NbDiracs_fin(3,k)+(min(index2(index2 > 0))==length(DiracsLocations)); 
  
   if nb == 1,  
    svd_ampl(:,k) = diag(S);
   else
    svd_ampl(:,k) = svd_ampl(:,k)+diag(S);
   end 
  end
 end
 svd_ampl = svd_ampl/Nbreal;
 NbDiracs_fin = NbDiracs_fin/Nbreal;
end
