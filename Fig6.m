%SNR = -10
close all;
par1 = 5:20;

[Non_Detect1] = test_detect([100 120],0,-10,par1,30);
[Non_Detect2] = test_detect([100 140],0,-10,par1,30);

plot(par1,Non_Detect1,par1,Non_Detect2,'--');
hold on;
%SNR = 0
[Non_Detect11] = test_detect([100 120],0,0,par1,30);
[Non_Detect21] = test_detect([100 140],0,0,par1,30);

plot(par1,Non_Detect11,'-d',par1,Non_Detect21,'->');
