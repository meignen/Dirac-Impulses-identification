%SNR = -10
par1 = 3:15;

[Non_Detect1] = test_detect([100 140 300],-5,-10,par1,30);
[Non_Detect2] = test_detect([100 140 300],-1,-10,par1,30);
[Non_Detect3] = test_detect([100 140 300],0,-10,par1,30);
[Non_Detect4] = test_detect([100 140 300],1,-10,par1,30);
[Non_Detect5] = test_detect([100 140 300],5,-10,par1,30);

plot(par1,Non_Detect1,par1,Non_Detect2,'--',par1,Non_Detect3,'-.',...
    par1,Non_Detect4,':',par1,Non_Detect5,'-s');

%SNR = 0
[Non_Detect1] = test_detect([100 140 300],-5,0,par1,30);
[Non_Detect2] = test_detect([100 140 300],-1,0,par1,30);
[Non_Detect3] = test_detect([100 140 300],0,0,par1,30);
[Non_Detect4] = test_detect([100 140 300],1,0,par1,30);
[Non_Detect5] = test_detect([100 140 300],5,0,par1,30);

figure
plot(par1,Non_Detect1,par1,Non_Detect2,'--',par1,Non_Detect3,'-.',...
    par1,Non_Detect4,':',par1,Non_Detect5,'-s');

%SNR = 10
[Non_Detect1] = test_detect([100 140 300],-5,10,par1,30);
[Non_Detect2] = test_detect([100 140 300],-1,10,par1,30);
[Non_Detect3] = test_detect([100 140 300],0,10,par1,30);
[Non_Detect4] = test_detect([100 140 300],1,10,par1,30);
[Non_Detect5] = test_detect([100 140 300],5,10,par1,30);

figure
plot(par1,Non_Detect1,par1,Non_Detect2,'--',par1,Non_Detect3,'-.',...
    par1,Non_Detect4,':',par1,Non_Detect5,'-s');
