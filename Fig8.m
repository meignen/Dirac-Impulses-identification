  sigma_w =40;   
  h = exp(-pi*(-1249:1250).^2/sigma_w^2);
  F = h';  
  Fc = toeplitz(F',[F(1) (F(2500:-1:2))']);
  
  LenF = length(F);
  hold on;
  for k = 0:7,
   DiracsLocations = [100 200-k*10];
   NbDiracs = length(DiracsLocations);
   DiracsWeights = ones(1,NbDiracs);
   s = zeros(LenF,1);
   s(DiracsLocations) = DiracsWeights;
  
   sf = real(ifft(fft(s).*fft(F)));
   
   %test MP
   mpdict = Fc;
   [yfit,r,coef,iopt,qual] = wmpalg('BMP',sf,mpdict,'itermax',NbDiracs);
   [X,I] = sort(iopt);
   
   %if one uses OMP, the results are identical
   %mpdict = Fc;
   %[yfit,r,coef,iopt,qual] = wmpalg('OMP',sf,mpdict,'itermax',NbDiracs);
   %[X1,I] = sort(iopt); 
   
   %% Non negative Matching Pursuit
   kk=1;
   r=sf;
   d=Fc'*r;
   F1=[];
   iopt=zeros(NbDiracs,1);
   while kk<=NbDiracs & max(d)>0
    [a b]=max(d);
    F1=[F1 Fc(:,b)];
    iopt(kk)=b;
    thet=lsqlin(F1,sf,[],[]);%,[],[],zeros(kk,1),[]);
    r=sf-F1*thet;
    d=Fc'*r;
    kk=kk+1;
   end
   coef=thet;
   
   [X2,I] = sort(iopt);
    
   plot(k,DiracsLocations(1),'o',k,DiracsLocations(2),'o',...
        k,X(1),'*',k,X(2),'*',k,X2(1),'s',k,X2(2),'s','LineWidth',4,'MarkerSize',20);
  end
  hold off;