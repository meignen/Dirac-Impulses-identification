function [FoundDiracsWeights,FoundDiracsWeights_locknown,FoundDiracsWeights_estime]...
                                             =test_weights(DiracsLocations,SNR,del,T)
  
 %% filter design
 
 sigma_w = 40;
 prec = 10^(-3);
 h = exp(-pi*(-1249:1250).^2/sigma_w^2);
 ffth = abs(fft(h));
 M0 = floor(sum((ffth/ffth(1) > prec))/2);
 F = h';  
 LenF = length(F);
 Fc = toeplitz(F',[F(1) (F(2500:-1:2))']);
 
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
 
 FoundDiracsWeights_locknown  = zeros(length(SNR),NbDiracs);
 FoundDiracsWeights_estime    = zeros(length(SNR),NbDiracs);
 FoundDiracsWeights           = zeros(length(SNR),NbDiracs);
  
 Nbreal = 200; 
 angle_TZroot = zeros(NbDiracs,1);
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
     
      while (abs(m_filt(k1)) <= 3.5*gamma_estime),
       X1 = [X1 k1];
       k1 = k1-1;
      end
      X = setdiff(X,X1);
      Xcur = setdiff(Xcur,X1);
      YY = Xcur(Xcur > 3.5*gamma_estime);
      if (length(YY) < 10)
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
    
   %test whether we have detected a the Dirac pulses
   n_detect = 0;
   for p = 1:NbDiracs,
    if abs(DiracLocations_estime(p)-DiracsLocations(p)) < 10    
     n_detect = n_detect+1;
    end 
   end
     
   %assume two Diracs
     
   %locations known
   a2i = Fc(:,[DiracsLocations(1) , DiracsLocations(2)])\m;
   FoundDiracsWeights_locknown(k,:) = FoundDiracsWeights_locknown(k,:)+...
                                        abs(a2i'-DiracsWeights);
  
   %error on weights using the estimated locations
   if n_detect == NbDiracs
    a2i = Fc(:,[DiracLocations_estime(1) , DiracLocations_estime(2)])\m;
    FoundDiracsWeights_estime(k,:) = FoundDiracsWeights_estime(k,:)+ abs(a2i'-DiracsWeights);
    nb_detect(k) = nb_detect(k)+1;
   end
     
   %minimisation of the error in the vicinity of the true locations, in
   %case of good detection
   if (n_detect == NbDiracs)
    [XX,YY] = meshgrid(DiracLocations_estime(1)+(-del:del),DiracLocations_estime(2)+(-del:del));
    error = inf;
    for p = 1:2*del+1,
     for q = 1:2*del+1, 
      if ((XX(p,q))<= LenF) && ((YY(p,q))<= LenF) 
       if (rank(Fc(:,[XX(p,q), YY(p,q)])) == 2)  
        a2i = Fc(:,[XX(p,q), YY(p,q)])\m;
        err = norm(m-Fc(:,[XX(p,q),YY(p,q)])*a2i);
        if err  < error
         error   = err;
         val_opt = a2i;
        end
       end
      end 
     end
    end
    FoundDiracsWeights(k,:) = FoundDiracsWeights(k,:)+ abs(val_opt'-DiracsWeights);
   end 
  end
 end
 FoundDiracsWeights_locknown  = FoundDiracsWeights_locknown/Nbreal;
 for k=1:length(SNR)
  FoundDiracsWeights_estime(k,:)  = FoundDiracsWeights_estime(k,:)./nb_detect(k);
  FoundDiracsWeights(k,:)         = FoundDiracsWeights(k,:)./nb_detect(k);
 end
end 