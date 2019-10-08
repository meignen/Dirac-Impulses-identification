
%Draws Fig 2 of the paper
close all 

%% Creation stream of Dirac
LenF = 2500;
NbDiracs = 2; 
DiracsWeights = [1 1];
DiracsLocations = [100 160];
s = zeros(LenF,1);
s(DiracsLocations) = DiracsWeights;
 
%% signal definition and frequency truncated version

sigma_w = 40;
h = exp(-pi*(-1249:1250).^2/sigma_w^2);
ffth = abs(fft(h));
M0 = floor(sum((ffth/ffth(1) > 10^(-3)))/2);
F = h';

sf = real(ifft(fft(s).*fft(F))); %F est supposé périodique, on fait convolution circulaire...
 
shat = 1/LenF*fft(sf);
sH = [conj(shat(M0+1:-1:2));shat(1:M0+1)];
shat_tronc = zeros(length(sf),1);
shat_tronc(1:M0+1)= shat(1:M0+1);
shat_tronc(end-M0+1:end)= shat(end-M0+1:end);
sf_filt = LenF*ifft(shat_tronc); %the reconstructed frequency truncated signal 

%% definition of phi and fequency truncation

phihat = 1/LenF*fft(F);
phiH   = zeros(length(phihat),1);
Dphi   = diag(phiH(1:2*M0+1));   

noise = randn(length(sf),1);
m     = sigmerge(sf,noise,-10);
mhat  = 1/LenF*fft(m);

mH = [conj(mhat(M0+1:-1:2));mhat(1:M0+1)];
   
mhat_tronc = zeros(length(mhat),1);
mhat_tronc(1:M0+1) = mhat(1:M0+1);
mhat_tronc(end-M0+1:end) = mhat(end-M0+1:end);
  
%% frequency truncation of m
m_filt = real(LenF*ifft(mhat_tronc));

%% estimation of the remaining noise level
gamma_estime = median(abs(m_filt-median(m_filt)))/0.6745;

%% filtering oscillations using piecewise cubic polynomial interpolation
X = (1:LenF).*((m_filt > 4*gamma_estime)'); 
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
m_filt1 = m_filt;
m_filt = pchip(abscissae,ordinates,1:LenF);

plot(X-LenF/2,m_filt(X)','*',Y-LenF/2,zeros(1,length(Y)),'o',1-LenF/2:LenF/2,m_filt,1-LenF/2:LenF/2,m_filt1,':',1-LenF/2:LenF/2,sf_filt,'--')

%% filtering oscillations using piecewise cubic polynomial interpolation
k = 1;
X = [];
m_filt = m_filt1;
while (k <= LenF)
 if (k<= LenF)
  while (abs(m_filt(k)) <= 4*gamma_estime),
   k = k+1;
   if (k > LenF)
    break;   
   end
  end
  
  if (k <= LenF)
   Xcur = [];
   while (m_filt(k) > 0.2*gamma_estime),
    X = [X k];
    Xcur = [Xcur k];
    k = k+1;
    if (k > LenF),
     break;
    end 
   end
  
   %we move backward to remove the points that were wrongly added
   if (k <= LenF)
    k1 = k;
    X1 = [];
    while (m_filt(k1) <= 4*gamma_estime),
     X1 = [X1 k1];
     k1 = k1-1;
    end
    X = setdiff(X,X1);
    Xcur = setdiff(Xcur,X1);
    YY = Xcur(Xcur > 4*gamma_estime);
    
    if length(YY) <= 10
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
 m_filt1 = m_filt;
 m_filt = pchip(abscissae,ordinates,1:LenF);

 figure 
 plot(X-LenF/2,m_filt(X)','*',Y-LenF/2,zeros(1,length(Y)),'o',1-LenF/2:LenF/2,m_filt,1-LenF/2:LenF/2,m_filt1,':',1-LenF/2:LenF/2,sf_filt,'--')
 