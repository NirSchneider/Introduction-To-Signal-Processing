close all;
clc;

disp('316098052, 315325654');
% Question 2

%% D

load('data2020.mat');
w=(0:1:length(y)-1);
Fs=44100;

x_test = [x_test, zeros(1,length(y_test)-length(x_test))];
X_test=fft(x_test,length(y));

Yz=fft(y_z);
Y_test=fft(y_test);

Hrec=(X_test)./(Y_test-Yz);
hrec=ifft(Hrec);

plot(w,abs(Hrec));
title('H_{rec}[k]'); xlabel('k'); ylabel('H_{rec}[k]');


figure;
plot(w,hrec);
title('h_{rec}[n]'); xlabel('n'); ylabel('h_{rec}[n]');

%% E
Y=fft(y);
Y0=Y-Yz;
Xrec=Hrec.*Y0;
xrec=ifft(Xrec);

figure;
plot(w,abs(Xrec));
title('X_{rec}[k]'); xlabel('k'); ylabel('X_{rec}[k]');

figure;
plot(w,xrec,[128241 128241],[-0.8 0.8],'--k');
title('x_{rec}[n]'); xlabel('n'); ylabel('x_{rec}[n]');

%% F

soundsc(xrec,Fs);
