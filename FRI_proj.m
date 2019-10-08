function [FoundDiracsLocation1,F] = FRI_proj(DiracsLocations,cas,cas1,sigma_w,prec,SNR,L,down_samp)
 
 %INPUT
 %DiracsLocations: the locations of Dirac pulses
 %cas            : if 0 no denoising algorithm, 
 %                 if 1 denoising algorithm (no downsampling)
 %                 if 2 denoising algorithm (with downsampling)
 %cas1           : if 0 the filter used is that of the LIDAR
 %                 if 1 the filter is the Gaussian filter
 %                 if 2 the filter is the periodic Dirichlet filter
 %OUTPUT
 %FoundDiracLocation1: the estimated Dirac pulses location
 %F                  : the filter F used
 
 if (cas1 == 0)
  folder = '/Users/Sylvain/Documents/home_kora/articles/Encour/Yoann_Sylvain/new_data_40m';
  fullMatFileName = fullfile(folder,  'F_R50_F500_325_100s_HiVis.mat');
  if ~exist(fullMatFileName, 'file')
   message = sprintf('%s does not exist', fullMatFileName);
   uiwait(warndlg(message));
  else
   load(fullMatFileName);
  end   
  F  = F(:,100);
  ffth = abs(fft(F));
  M0 = floor(sum((ffth/ffth(1) > prec))/2);
  
 elseif (cas1 == 1) 
  h = exp(-pi*(-1249:1250).^2/sigma_w^2);
  ffth = abs(fft(h));
%   X=abs(fftshift(fft(h)));
%   X(1250)
%   X(1280)
%   plot(-1249:1250,abs(fftshift(fft(h))));
%   pause
  M0 = floor(sum((ffth/ffth(1) > prec))/2);
  F = h';  
 end
 
 if cas1 < 2;
  LenF = length(F);
 else
  LenF =2500;
  M0 = 92;
 end
 
 %% Creation stream of Diracs
 
 NbDiracs = length(DiracsLocations);
 DiracsWeights = ones(1,NbDiracs);
 s = zeros(LenF,1);
 s(DiracsLocations) = DiracsWeights;
 
 %% signal definition

 if cas1 < 2
  sf = real(ifft(fft(s).*fft(F)));   
 else
  %periodic Dirichlet kernel
  fftF = zeros(LenF,1); 
  fftF(1:M0+1) = ones(M0+1,1);
  fftF(end:-1:end-M0+1) = ones(M0,1);
  sf = real(ifft(fft(s).*fftF));
  F = 0;
 end 
 
 %% definition of phi and fequency truncation

 if cas1 < 2,
  phihat = 1/LenF*fft(F);
  phiH  = [conj(phihat(M0+1:-1:2));phihat(1:M0+1)];  
 else
  phiH  = ones(2*M0+1,1);
 end 
 
 FoundDiracsLocation1 = zeros(length(SNR),NbDiracs);
 
 Nbreal = 200; 
 angle_TZroot = zeros(NbDiracs,1);
 DiracLocations_estime = zeros(NbDiracs,down_samp);
 add_twopi = (0:down_samp-1)*2*pi;
 F_dec = zeros(LenF,1);
 res_detect = zeros(NbDiracs,down_samp*NbDiracs);
 val_Dir = zeros(1,NbDiracs);
 
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
   
   if (cas == 2)||(cas == 1)
    %estimation of the noise level
    gamma_estime = median(abs(m_filt-median(m_filt)))/0.6745;
  
    %filtering oscillations using piecewise cubic polynomial interpolation
    X = (1:LenF).*((abs(m_filt) > 5*gamma_estime)'); 
    X = X(X>0);
    Y = (1:LenF).*((abs(m_filt) < 0.2*gamma_estime)');
    Y = Y(Y>0);
    abscissae = [X Y];
    if ismember(LenF,abscissae)
     ordinates =  [m_filt(X)' zeros(1,length(Y))];
    else
     abscissae = [X Y LenF];
     ordinates =  [m_filt(X)' zeros(1,length(Y)) 0];
    end
    
    %piecewise cubic interpolation
    m_filt = pchip(abscissae,ordinates,1:LenF);
    m_filt = m_filt(:);
    mhat = 1/LenF*fft(m_filt);
  
    mH = [conj(mhat(M0+1:-1:2));mhat(1:M0+1)];
   end 
   
   y  = mH./phiH;
   
   %we select the indices 
   indice = zeros(2*L+1,1);
   indice(L+1:2*L+1) = M0+1:down_samp:M0+1+L*down_samp;
   indice(1:L) = M0+1-L*down_samp:down_samp:M0+1-down_samp;
  
   y_t = y(indice);
   ind = NbDiracs+1;
   A = toeplitz(y_t(ind:ind+2*L-NbDiracs),y_t(ind:-1:ind-NbDiracs)); 
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
   
   if (cas == 0)||(cas == 1)
    DiracLocations = sort(round(angle_TZroot/(2*pi)*LenF)+1);
    if nb == 1
     FoundDiracsLocation1(k,:) = DiracLocations;
    else     
     X = zeros(NbDiracs,1);
     X(:) = FoundDiracsLocation1(k,:);
     FoundDiracsLocation1(k,:) = X+ DiracLocations;
    end
   end
    
   if cas == 2,
    Ang_test = zeros(NbDiracs,down_samp); 
    for p = 1:NbDiracs,    
     Ang_test(p,:) = angle_TZroot(p)*ones(1,down_samp)+add_twopi;        
     DiracLocations_estime(p,:) = round(Ang_test(p,:)/(2*pi*down_samp)*LenF)+1;
    end
    DiracsEstime = DiracLocations_estime(:);
   
    for p = 1:NbDiracs,
     for q = 1:down_samp*NbDiracs,
      F_dec(1:DiracsEstime(q)-1) = F(end-DiracsEstime(q)+2:end);
      F_dec(DiracsEstime(q):end) = F(1:end-DiracsEstime(q)+1); 
      res_detect(p,q) = sum(m.*F_dec);
     end
     [val,ind_NbDiracs] = sort(res_detect(p,:),'descend');
     if p == 1,
      val_Dir(p) = DiracsEstime(ind_NbDiracs(1));
     else
      ii = 1;   
      while  min(abs(val_Dir(p-1)-DiracsEstime(ind_NbDiracs(ii)))) < 10
       ii = ii +1;
      end 
      val_Dir(p) = DiracsEstime(ind_NbDiracs(ii));
     end
    end
    
    if nb == 1
     FoundDiracsLocation1(k,:) = sort(val_Dir);
    else
     FoundDiracsLocation1(k,:) = FoundDiracsLocation1(k,:)+sort(val_Dir);
    end
   end 
  end
 end 
 FoundDiracsLocation1 = FoundDiracsLocation1/Nbreal;
end 
