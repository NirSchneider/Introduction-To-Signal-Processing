close all;
clc;

disp('316098052, 315325654');
% Question 3

%% B
T=10;
cn=-20:20;
t=(0:1/100:T);
fi=exp((1j*2*pi*t.'*cn)/T);  
n=(0:19)';
xsi=(rectangularPulse((t./(T/20))-n-0.5))';

f=(4*cos((4*pi*t)/T)+sin((10*pi*t)/T))';
g=(2*sign(sin((6*pi*t)/T))-4*sign(sin((4*pi*t)/T)))';

Cn_f_fi=coff(f,fi,T);
Cn_f_xsi=coff(f,xsi,T);
Cn_g_fi=coff(g,fi,T);
Cn_g_xsi=coff(g,xsi,T);

%% C
SUM_f_fi_temp=fi*Cn_f_fi;
SUM_f_xsi=xsi*Cn_f_xsi;
SUM_g_fi_temp=fi*Cn_g_fi;
SUM_g_xsi=xsi*Cn_g_xsi;

for i=1:1001
    SUM_f_fi(i)=SUM_f_fi_temp(1002-i);
end
for i=1:1001
    SUM_g_fi(i)=SUM_g_fi_temp(1002-i);
end

figure;
plot(t,f,'b',t,SUM_f_fi,'m');
legend('f','f_r_e_c'); title('Reconstruction f(t)  fron  fi_n(t)'); ylabel('f(t)'); xlabel('t [sec]');

figure;
plot(t,f,'b',t,SUM_f_xsi,'g');
legend('f','f_r_e_c'); title('Reconstruction f(t)  fron  xsi_n(t)'); ylabel('f(t)'); xlabel('t [sec]');

figure;
plot(t,g,'b',t,SUM_g_fi,'m');
legend('g','g_r_e_c'); title('Reconstruction g(t)  fron  fi_n(t)'); ylabel('g(t)'); xlabel('t [sec]');

figure;
plot(t,g,'b',t,SUM_g_xsi,'g');
legend('g','g_r_e_c'); title('Reconstruction g(t)  fron  xsi_n(t)'); ylabel('g(t)'); xlabel('t [sec]');

%% functions
function C = coff(vec,Mat,T)
    [m,n]=size(Mat);
    t=linspace(0,T,m);
    complex=conj(Mat);
    
    for i=1:n
       
        M(i)=(trapz(t,(vec.').*(complex(:,i).'))) ./ (trapz(t,(Mat(:,i).').*(complex(:,i).')));
    end
    C=M';
end