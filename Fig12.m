close all;
sigma_w = 40;
prec = 10^(-3);
h_orig = exp(-pi*(-1249:1250).^2/sigma_w^2);

gaussian = @(x) (1/sqrt((2*pi))*exp(-x.^2/2));
skewedgaussian = ...
 @(x,alpha) 2*sqrt(2*pi)*gaussian(sqrt(2*pi)/sigma_w*x).*...
                  normcdf(alpha*sqrt(2*pi)/sigma_w*x);

h1 = skewedgaussian(-1249:1250,1);
h5 = skewedgaussian(-1249:1250,5);
h_1 = skewedgaussian(-1249:1250,-1);
h_5 = skewedgaussian(-1249:1250,-5);
plot(1150:1350,h_orig(1150:1350),1150:1350,h1(1150:1350),'--',1150:1350,h5(1150:1350),'-.',...
               1150:1350,h_1(1150:1350),'-<',1150:1350,h_5(1150:1350),'->')