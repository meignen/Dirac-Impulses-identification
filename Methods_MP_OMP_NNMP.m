function [FoundDiracsLocations_MP,FoundDiracsLocations_NNMP,...
          FoundDiracsWeights_MP,FoundDiracsWeights_NNMP]...
           = Methods_MP_OMP_NNMP(DiracsLocations,cas,SNR)

 %% Lecture filtre phi
 if (cas == 0)
  folder = '/Users/Sylvain/Documents/home_kora/articles/Encour/Yoann_Sylvain/new_data_40m';
  fullMatFileName = fullfile(folder,'F_R50_F500_325_100s_HiVis.mat');
  if ~exist(fullMatFileName, 'file')
   message = sprintf('%s does not exist', fullMatFileName);
   uiwait(warndlg(message));
  else
   load(fullMatFileName);
  end
  Fc = F(:,101:2200);
  F  = F(:,100);
 else
  sigma_w =40;   
  h = exp(-pi*(-1249:1250).^2/sigma_w^2);
  F = h';  
  Fc = toeplitz(F',[F(1) (F(2500:-1:2))']);
 end
 
 %% Creation stream of Dirac
 LenF = length(F);
 
 NbDiracs = length(DiracsLocations);
 DiracsWeights = ones(1,NbDiracs);
 s = zeros(LenF,1);
 s(DiracsLocations) = DiracsWeights;
 

 %% Filtrage par Phi
 sf = real(ifft(fft(s).*fft(F)));
 %sf=Fc(:,DiracsLocations)*DiracsWeights';

 %% Add noise
 Nbreal = 50;

 FoundDiracsLocations_MP  = zeros(length(SNR),NbDiracs);
 %FoundDiracsLocations_OMP = zeros(NbDiracs,length(SNR));
 FoundDiracsLocations_NNMP = zeros(length(SNR),NbDiracs);
 FoundDiracsWeights_MP   = zeros(length(SNR),NbDiracs);
 %FoundDiracsWeights_OMP  = zeros(NbDiracs,length(SNR));
 FoundDiracsWeights_NNMP = zeros(length(SNR),NbDiracs);

for nb =1:Nbreal
  for k= 1:length(SNR)  
   noise = randn(length(sf),1);   
   m = sigmerge(sf,noise,SNR(k));

   %% Mathcing pursuit
   mpdict = Fc;
   [yfit,r,coef,iopt,qual] = wmpalg('BMP',sf,mpdict,'itermax',NbDiracs);
    
   [X,I] = sort(iopt);   
   X
   pause
   FoundDiracsLocations_MP(k,:) = FoundDiracsLocations_MP(k,:)+abs(X -DiracsLocations);
   FoundDiracsWeights_MP(k,:)   = FoundDiracsWeights_MP(k,:)+ abs((coef(I))'-DiracsWeights); 
  
   %% Orthogonal Matching pursuit
%    mpdict = Fc;
%    [yfit,r,coef,iopt,qual] = wmpalg('OMP',m,mpdict,'itermax',NbDiracs);
%    [X,I] = sort(iopt);   
%    FoundDiracsLocations_OMP(:,k) = FoundDiracsLocations_OMP(:,k)+abs(X-DiracsLocations)';
%    FoundDiracsWeights_OMP(:,k)   = FoundDiracsWeights_OMP(:,k)+ abs(coef(I)-DiracsWeights');
   
   %% Non negative Matching Pursuit
   kk=1;
   r=m;
   d=Fc'*r;
   F1=[];
   iopt=zeros(NbDiracs,1);
   while kk<=NbDiracs & max(d)>0
    [a b]=max(d);
    F1=[F1 Fc(:,b)];
    iopt(kk)=b;
    thet=lsqlin(F1,m,[],[]);%,[],[],zeros(kk,1),[]);
    r=m-F1*thet;
    d=Fc'*r;
    kk=kk+1;
   end
   coef=thet;
   
   [X,I] = sort(iopt);   
   FoundDiracsLocations_NNMP(k,:) = FoundDiracsLocations_NNMP(k,:)+ abs(X'-DiracsLocations);
   FoundDiracsWeights_NNMP(k,:) = FoundDiracsWeights_NNMP(k,:)+ abs((coef(I))'-DiracsWeights);
 end
end

FoundDiracsLocations_MP = FoundDiracsLocations_MP/Nbreal;
%FoundDiracsLocations_OMP = FoundDiracsLocations_OMP/Nbreal;
FoundDiracsLocations_NNMP = FoundDiracsLocations_NNMP/Nbreal;

FoundDiracsWeights_MP = FoundDiracsWeights_MP/Nbreal;
%FoundDiracsWeights_OMP = FoundDiracsWeights_OMP/Nbreal;
FoundDiracsWeights_NNMP = FoundDiracsWeights_NNMP/Nbreal;
