close all;
clc;

disp('316098052, 315325654');
% Question 1

%% A

wm=3*pi;
t=0.2:1/100:3; %continuous time vector
x=4./(wm*pi*t.^2).*(sin(wm*t)).^2.*(cos(wm*t)).*(sin(2*wm*t)); %continuous signal

figure; %new figure window
plot(t,x,'k','LineWidth',2); %draw the continuous signal graph - plot(t,x,'LineStyle','-','Color','k','LineWidth',2);
title('Q 1,a - Signal x(t) - time domain'); xlabel('t [sec]','FontSize',12); ylabel('x(t)[V]','FontSize',12);

figure; %absolute value
plot(t,abs(x),'k','LineWidth',2); %draw the continuous signal graph - plot(t,x,'LineStyle','-','Color','k','LineWidth',2);
title('Q 1,a '); xlabel('t [sec]','FontSize',12); ylabel('x(t)[V]','FontSize',12);

%% B

w=(-17*pi:0.01:17*pi);
xf=1/1j.*(triangularPulse((w-wm)/(2*wm))-triangularPulse((w+wm)/(2*wm))+triangularPulse((w-3*wm)/(2*wm))-triangularPulse((w+3*wm)/(2*wm)));
plot(w,abs(xf));
title('Q 1,b - X(jw) magnitude '); xlabel('w [rad/sec]','FontSize',12); ylabel('| X(jw) |','FontSize',12);

%% C

ws=10*wm;
Ts=(2*pi)/ws;
ts=0.2:Ts:3;
xZOH=4./(wm*pi*ts.^2).*(sin(wm*ts)).^2.*(cos(wm*ts)).*(sin(2*wm*ts));
figure;
stairs(ts, xZOH, 'b');
hold on;
plot(t,x,'k','LineWidth',1);
plot(ts, xZOH,'r*');
legend('x(nTs)','x(t)','samples')
hold off;
title('Q1.c - x(nTs) ZOH - time domain');
ylabel('x(nTs)');
xlabel('t [sec]');

%% D

xp=xf+1/1j.*(triangularPulse((w-wm+ws)/(2*wm))-triangularPulse((w+wm+ws)/(2*wm))+triangularPulse((w-3*wm+ws)/(2*wm))-triangularPulse((w+3*wm+ws)/(2*wm)))+1/1j.*(triangularPulse((w-wm-ws)/(2*wm))-triangularPulse((w+wm-ws)/(2*wm))+triangularPulse((w-3*wm-ws)/(2*wm))-triangularPulse((w+3*wm-ws)/(2*wm)));
%2 sampling is enough because after those sampling its not the range of W

figure;
xfZOH=exp(-(1j*Ts/2)*w).*sinc(w/ws).*xp;
plot(w,abs(xfZOH));
title('Q1.d - |Xzoh(jw)| - frequency domain');
ylabel('Xzoh(jw)');
xlabel('w [rad]');

%% E

size=length(w);
H=zeros(1,size);

for s=1:size
    if (w(s)>-ws/2&&w(s)<ws/2)
        H(s)=exp((1j*pi*w(s))/ws)/(sinc(w(s)/ws));
    else
        H(s)=0;
    end
end

xfrec=H.*xfZOH;

xrec=zeros(1,281);

for k=1:281
xrec(k)=(1/(2*pi))*trapz(w,xfrec.*(exp(1j*w*t(k))));
end

figure;
plot(t,x);
hold on;
plot(t,xrec,'--');
hold off;
legend('x(t)','xrec(t)')
title('Q1.e - the recoverd signal');
ylabel('xrec(t)');
xlabel('t [sec]');

%% F

ws_new=9*wm;
Ts_new=(2*pi)/ws_new;
 Xzoh_F1=0;
for k=-2:2
   Xzoh_F1=Xzoh_F1+ 1/1i.*(triangularPulse((w+k*ws_new)/(6*pi)-1.5)-triangularPulse((w+k*ws_new)/(6*pi)+0.5)+triangularPulse((w+k*ws_new)/(6*pi)-0.5)-triangularPulse((w+k*ws_new)/(6*pi)+1.5));
end
Xzoh_F1=Xzoh_F1.*exp(-1i.*w.*Ts/2).*sinc(w/(2*pi/Ts_new));
Xrec_F1=Xzoh_F1./(sinc(w/ws_new)).*exp(1i.*w.*Ts_new/2).*(heaviside(w+(ws_new/2))-heaviside(w-(ws_new/2)));
xrec1=zeros(1,length(t));
for i=1:length(t)
    xrec1(i)=1/(2*pi)*trapz(w,Xrec_F1.*exp(1i*w*(t(i))));
end
plot(t,x,'b','LineWidth',2);
hold on;
plot(t,xrec1,'r','LineWidth',2);
legend('x(t)','xrec(t)-ws=9wm')
title('x(t) , xrec(t)');
ylabel('xrec(t)');
xlabel('t [sec]');
