close all;
clc;

disp('316098052, 315325654');
% Question 1

%% A

teta1=pi/10.25;
teta2=2*pi/5;
N=30;
n=0:(N-1);
s_n=2*cos(teta1*n);
v_n=3*sin(teta2*n);
x_n=v_n+s_n;

plot(n,x_n,'r-*',n,v_n,'b-o',n,s_n,'g-+');grid on;
title('The signals x[n],v[n],s[n]'); ylabel('x[n] s[n] v[n]'); xlabel('n'); legend('x[n]', 'v[n]', 's[n]');

% n_k = (2*pi/N)*n;
s_k=fft(s_n);
v_k=fft(v_n);
x_k=fft(x_n);

figure;
hold on;
stem(n,abs(x_k),'r-*');
stem(n,abs(v_k),'b-o');
stem(n,abs(s_k),'g-+');
hold off;
title('DFT of signals x[n],v[n],s[n]'); ylabel('X[k] V[k] S[k]'); xlabel('k'); legend('X[k]' ,'V[k]', 'S[k]');

%% B

xz_n=[x_n ,zeros(1,15)];
xz_k=fft(xz_n);
Nz=45;
nz=(0:N/Nz:(N-1/2));

figure;
hold on;
stem(n,abs(x_k),'r-*');
stem(nz,abs(xz_k),'b-o');
hold off;
title('X[k] , X_z[k_z]'); ylabel('X[k] , X_z[k_z]'); xlabel('k'); legend('X[k]' ,'X_z[k_z]');

%% C

N2=45;
nc=0:(N2-1);
s_n2=2*cos(teta1*nc);
v_n2=3*sin(teta2*nc);
x_n2=v_n2+s_n2;
x_k2=fft(x_n2);

figure
hold on;
stem(n,abs(x_k),'r-*');
stem(nz,abs(x_k2),'k-o');
% stem(nz,abs(xz_k),'g-o');
hold off;
title('X[k] , X_2[k_2]'); ylabel('X[k] , X_2[k_2]'); xlabel('k'); legend('X[k]' ,'X_2[k_2]');
% title('X[k] , X_2[k_2], X_z[k]'); ylabel('X[k] , X_2[k_2]'); xlabel('k'); legend('X[k]' ,'X_2[k_2], X_z[k]');

%% D

P_n=x_n*x_n';
P_k=x_k*x_k'/N;

P_zn=xz_n*xz_n';
P_zk=xz_k*xz_k'/N2;

%% E

h1=[1/3 ,1/3 ,1/3 zeros(1,29)];
H1=fft(h1);
ny=(0:1:31);

xzz_n=[x_n ,zeros(1,2)];
Xzz_k=fft(xzz_n);

figure;
hold on;
plot(ny,abs(H1)/max(abs(H1)),'r');
stem(ny,abs(Xzz_k)/max(abs(Xzz_k)),'b-o');
hold off;
title('X[k] and H_1(\theta)'); ylabel('X[k], H_1(\theta)'); xlabel('k'); legend('H_1(\theta)', 'X[k]');

Y=Xzz_k.*H1;
y=ifft(Y);

figure;
hold on;
stem(ny,abs(Y),'r-*', 'linewidth', 1);
stem(ny,abs(Xzz_k),'b-o');
hold off;
title('Y[k] , X[k]'); ylabel('Y[k] , X[k]'); xlabel('k'); legend('Y[k]' ,'X[k]');

figure;
hold on;
plot(n,y(1:30),'k-o', 'linewidth', 1);
plot(n,x_n,'r-*');
plot(n,s_n,'b-*');
plot(n,v_n,'g-+');
hold off;
title('y[n],x[n],s[n],v[n]'); ylabel('y[n],x[n],s[n],v[n]'); xlabel('n'); legend('y[n]','x[n]','s[n]','v[n]');

%% F

h2=[1 ,1 zeros(1,29)];
H2=fft(h2);

xzz2_n=[x_n ,0];
Xzz2_k=fft(xzz2_n);

ny2=(0:1:30);
figure;
hold on;
plot(ny2,abs(H2)/max(abs(H2)),'r');
stem(ny2,abs(Xzz2_k)/max(abs(Xzz2_k)),'b-o');
hold off;
title('H_2(\theta) and X[k]'); ylabel('H_2(\theta) , X[k]'); xlabel('k'); legend('H_2(\theta)','X[k]');
 
Y2=Xzz2_k.*H2;
y2=ifft(Y2);

figure;
hold on;
stem(ny2,abs(Xzz2_k),'b-o');
stem(ny2,abs(Y2),'r-*', 'linewidth', 1);
hold off;
title('Y_2[k] and X[k]'); ylabel('Y_2[k], X[k]'); xlabel('k'); legend('X[k]','Y_2[k]');

figure;
hold on;
plot(n,y2(1:30),'k-o', 'linewidth', 1);
plot(n,x_n,'r-*');
plot(n,s_n,'b-*');
plot(n,v_n,'g-+');
hold off;
title('y_2[n],x[n],s[n],v[n]'); ylabel('y_2[n],x[n],s[n],v[n]'); xlabel('n'); legend('y_2[n]','x[n]','s[n]','v[n]');

%% G
%% case 1
teta11=pi/2;
teta21=(8*pi)/5;

s_n1=2*cos(teta11*n);
v_n1=3*sin(teta21*n);
x_n1=v_n1+s_n1;


h1=[1/3 ,1/3, 1/3];
y11=conv(h1,x_n1);
h2=[1 ,1];
y12=conv(h2,x_n1);

figure
hold on;
plot(n,s_n1,'k-*', 'linewidth', 1);
plot(n,y11(1:30),'b-o');
plot(n,y12(1:30),'r-o');
hold off;
title('s[n],y_1[n],y_2[n] - \theta_1=\pi/2 \theta_2=(8*\pi)/5'); ylabel('s[n],y_1[n],y_2[n]'); xlabel('n'); legend('s[n]','y_1[n]','y_2[n]');
%% case 2
teta12=pi/2;
teta22=pi;

s_n2=2*cos(teta12*n);
v_n2=3*sin(teta22*n);
x_n2=v_n2+s_n2;

h1=[1/3 ,1/3, 1/3];
y12=conv(h1,x_n2);
h2=[1 ,1];
y22=conv(h2,x_n2);

figure
hold on;
plot(n,s_n2,'k-*', 'linewidth', 1);
plot(n,y12(1:30),'b-o');
plot(n,y22(1:30),'r-o');
hold off;
title('s[n],y_1[n],y_2[n] - \theta_-1=pi/2 \theta_2=\pi'); ylabel('s[n],y_1[n],y_2[n]'); xlabel('n'); legend('s[n]','y_1[n]','y_2[n]');


