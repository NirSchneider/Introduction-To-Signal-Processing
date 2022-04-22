close all;
clc;

disp('316098052, 315325654');
% Question 2

%% A
wa=7*pi;
wb=4*pi;
Ta = 2; %period time of x
t=0:0.01:Ta;
ts=0:2/15:28/15;

x=5*cos(wa*t)-3*sin(wb*t);
xs=5*cos(wa*ts)-3*sin(wb*ts);

figure
hold on;
plot(t,x);
plot(ts,xs,'r*');
legend('x(t)','x_s(t)'); 
ylabel('x(t), x_sample(t)'); 
xlabel('t [sec]'); 
title('uniform sampled signal & the original signal');
hold off;

%% B
k=1:15;
Fs=exp(1i*(ts.').*(k-8)*pi);
F=exp(1i*(t.').*(k-8)*pi);
a=inv(Fs)*xs.';
disp(a);

%% C
xrec=F*a;
figure;
plot(t,x,'m',t,xrec,'--g');
grid on;
legend('x(t)','x rec(t)'); 
title('The recoverd signal - uniform sampling'); 
ylabel('x(t),x rec(t)'); 
xlabel('t [sec]');
hold off;

%% D
tsrand=Ta*rand(1,15);
xsrand=5*cos(wa*tsrand)-3*sin(wb*tsrand);
figure;
hold on;
plot(t,x);
plot(tsrand,xsrand,'r*');
legend('x(t)','x_srand(t)'); 
title('not-uniform sample signal & the original signal'); 
ylabel('x(t), x_srand(t)'); 
xlabel('t [sec]');
hold off;

Fsrand=exp(1i*(tsrand.').*(k-8)*pi);
arand=inv(Fsrand)*xsrand.';
xrecrand=F*arand;
figure;
plot(t,x,'m',t,xrecrand,'--g');
legend('x(t)','x-rec-rand(t)'); 
title('The recoverd signal - not-uniform sampling'); 
ylabel('x(t), x-rec-rand(t)'); 
xlabel('t [sec]');
hold off;

%% E
ts1 = zeros(1,15);
tsrand1 = zeros(1,15);
for l=1:15
    ts1(l)=ts(l)+0.01*rand(1);
    tsrand1(l)=tsrand(l)+0.01*rand(1);
end

Fs1=exp(1i*((ts1).').*(k-8)*pi);
a1=inv(Fs1)*xs.';
xrec1=F*a1;
figure;
plot(t,x,'m',t,xrec1,'--g',ts,xs,'r*');
grid on;
legend('x(t)','xrec(t)','x_s(t)');
title('Reconstruction from uniform sampling with random error'); 
ylabel('x(t)'); 
xlabel('t [sec]');
hold off;

Fsrand1=exp(1i*((tsrand1).').*(k-8)*pi);
arand1=inv(Fsrand1)*xsrand.';
xrecrand1=F*arand1;
figure;
plot(t,x,'m',t,xrecrand1,'--g',tsrand,xsrand,'r*');
grid on;
legend('x(t)','xrecrand1','x_srand(t)'); 
title('Reconstruction from not-uniform sampling with random error'); 
ylabel('x(t)');
xlabel('t [sec]');
hold off;

C1=cond(Fs);
C2=cond(Fsrand);
C3=cond(Fs1);
C4=cond(Fsrand1);

%% F
ts2=(0:1/20:39/20);
tsrand2=2*rand(1,40);
Fs2=exp(1i*((ts2).').*(k-8)*pi);
Fsrand2=exp(1i*((tsrand2).').*(k-8)*pi);

ts2_err = zeros(1,40);
tsrand2_err = zeros(1,40);
for l=1:40
    ts2_err(l)=ts2(l)+0.01*rand(1);
    tsrand2_err(l)=tsrand2(l)+0.01*rand(1);
end

%uniform Sampling
xs2=5*cos(wa*ts2)-3*sin(wb*ts2);
Fs2_err=exp(1i*((ts2_err).').*(k-8)*pi);
a2=inv((Fs2_err')*Fs2_err)*(Fs2_err')*xs2.';
xrec2=F*a2;
figure;
plot(t,x,'m',t,xrec2,'--g',ts2,xs2,'r*');
grid on;
legend('x(t)','xrec(t)','x_s(t)');
title('Reconstruction from uniform sampling with random error - 40 sampling points'); 
ylabel('x(t)'); 
xlabel('t [sec]');
hold off;


%non-uniform Sampling
xsrand2=5*cos(wa*tsrand2)-3*sin(wb*tsrand2);
Fsrand2_err=exp(1i*((tsrand2_err).').*(k-8)*pi);
arand2=inv((Fsrand2_err')*Fsrand2_err)*(Fsrand2_err')*xsrand2.';
xrecrand2=F*arand2;
figure;
plot(t,x,'m',t,xrecrand2,'--g',tsrand2,xsrand2,'r*');
grid on;
legend('x(t)','xrecrand1','x_srand(t)'); 
title('Reconstruction from non-uniform sampling with random error - 40 sampling points'); 
ylabel('x(t)');
xlabel('t [sec]');
hold off;

C5=cond(Fs2);
C6=cond(Fsrand2);
C7=cond(Fs2_err);
C8=cond(Fsrand2_err);