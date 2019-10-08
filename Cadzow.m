function Cadzow(DiracsLocations,sigma_w,prec,SNR,L)
  
  %creation of the filter
  h = exp(-pi*(-1249:1250).^2/sigma_w^2);
  ffth = abs(fft(h));

  M0 = floor(sum((ffth/ffth(1) > prec))/2);
  F = h';
  LenF = length(F);
   
  %creation of the stream of Diracs
  NbDiracs = length(DiracsLocations);
  DiracsWeights = ones(1,NbDiracs);
  s = zeros(LenF,1);
  s(DiracsLocations) = DiracsWeights;	
 
  %definition of the filtered signal
  sf = real(ifft(fft(s).*fft(F))); 
  
  %the truncated 
  phihat = 1/LenF*fft(F);
  phiH  = [conj(phihat(M0+1:-1:2));phihat(1:M0+1)];  
  
  noise = randn(length(sf),1);
  m = sigmerge(sf,noise,SNR);
  mhat = 1/LenF*fft(m);
   
  %the truncated values of mH between frequencies -M0 and M0
  mH = [conj(mhat(M0+1:-1:2));mhat(1:M0+1)];
  
  y  = mH./phiH; %the signal y with left 2*M0+1 
   
   %we select the central indices 
   indice = zeros(2*L+1,1);
   indice(L+1:2*L+1) = M0+1:M0+1+L;
   indice(1:L)       = M0+1-L:M0;
  
   K     = length(DiracsLocations);
   vhatm = y(indice);
   %ind = NbDiracs+1;
   %A = toeplitz(y_t(ind:ind+2*L-NbDiracs),y_t(ind:-1:ind-NbDiracs)); 
   
   %[U,S,W] = svd(A);
   
   N = 2*L+1;  %we consider a square matrix
   P = L;				%the matrices have size N-P x P+1. K<=P<=M required.
   Nbiter=50;			%number of iterations. 
   Tdenoised=toeplitz(vhatm(P+1:N),vhatm(P+1:-1:1)); %the noisy matrix is the initial estimate
   
   for iter=1:Nbiter
	[U S V]=svd(Tdenoised);
	Tdenoised=U(:,1:K)*S(1:K,1:K)*(V(:,1:K))';	%SVD truncation -> Tdenoised has rank K
	Tdenoised=Toeplitzation(Tdenoised);			% -> Tdenoised is Toeplitz	
   end
   
   %denoised vector
   y_denoised = zeros(2*L+1,1);
   y_denoised(L+1:2*L+1) = Tdenoised(:,1);
   y_denoised(1:L)       = Tdenoised(1,end:-1:2);
   
   %comparaison avec ytrue
   mhat_true = 1/LenF*fft(sf);
   
   %the truncated values of mH between frequencies -M0 and M0
   mHtrue = [conj(mhat_true(M0+1:-1:2));mhat_true(1:M0+1)];
   ytrue  = mHtrue./phiH;
   plot(1:2*L+1,real(ytrue(M0+1-L:M0+1+L)),1:2*L+1,real(y_denoised),'--');
   pause
end

function Matres=Toeplitzation(Mat)
   %this function returns a Toeplitz matrix, closest to Mat for the Frobenius norm.
   %this is done by simply averaging along the diagonals.
   [height,width]=size(Mat);  %height>=width required
    Matres=Mat;
	for indexcol=2:width
		valdiag=0;
		valdiag2=0;
		for indexrow=1:width-indexcol+1
		 valdiag=valdiag+Mat(indexrow,indexcol-1+indexrow);
		 valdiag2=valdiag2+Mat(height-indexrow+1,width-indexcol+2-indexrow);
		end
		valdiag=valdiag/(width-indexcol+1);
		valdiag2=valdiag2/(width-indexcol+1);
		for indexrow=1:width-indexcol+1
			Matres(indexrow,indexcol-1+indexrow)=valdiag;
			Matres(height-indexrow+1,width-indexcol+2-indexrow)=valdiag2;
		end
	end
	for indexcol=1:height-width+1
		valdiag=0;
		for indexrow=1:width
			valdiag=valdiag+Mat(indexcol+indexrow-1,indexrow);
		end
		valdiag=valdiag/width;
		for indexrow=1:width
			Matres(indexcol+indexrow-1,indexrow)=valdiag;
		end
    end
  end		%of the function	


