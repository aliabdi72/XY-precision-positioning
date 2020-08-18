% In the name of Allah
% Ali Abdi
% 1397-09-22
% Thesis
clear all
close all
clc
% Main Code
%% ------------------------------------- Input
%Input Desired Position
Xd=input('please enter Desired X (um) :   ')*1e-6
Yd=input('please enter Desired Y (um):   ')*1e-6
%% ------------------------------------- Parameters
% In common  
Ap=2.5e-5;
Ep=2.216e10;
ROp=8000;
lp=2e-2;
Mp=ROp*Ap*lp;
g=9.81;
ALFp=4.297e-1;
BETp=3.438e-2;
GAMp=-2.865e-3;
de=1.24e-7;
cp=25;
kp=Ep*Ap/lp;
n=1;
vs=0.001;
ALF=2e10;
BET=0.2e10;
GAM=0.3e10;
DEL=2;
SIG1=100;
SIG2=0.4;
% X
t0x=0e-3;
t1x=0.5e-3;
t2x=2.5e-3;
t3x=3e-3;
t4x=8e-3;       % Rest Time
dtx=1e-7;
ssx1=(2*t1x-t3x)/(t3x-t0x);
ssx2=(2*t2x-t3x)/(t3x-t0x);
T1x=[-t3x/2:dtx:+t3x/2];
datx1=length(T1x);
Tx=[0:dtx:t4x];
datx=length(Tx);
Vmaxxf=54;
Vmaxxb=49;
Volxf=zeros(1,datx);  %Forward
Volxb=zeros(1,datx);  %Backward
Volxf(1:datx1)=Vmaxxf*tripuls(T1x,t3x,ssx1);  %Forward
Volxb(1:datx1)=Vmaxxb*tripuls(T1x,t3x,ssx2);  %Backward
m1x=400e-3;
m2x=300e-3; 
MUkx=0.5;
MUsx=0.6;
Fnx=(m1x+m2x+Mp)*g;
Fsx=MUsx*Fnx;                      
Fcx=MUkx*Fnx;
% Y
t0y=0e-3;
t1y=0.5e-3;
t2y=2.5e-3;
t3y=3e-3;
t4y=8e-3;       % Rest Time
dty=1e-7;
ssy1=(2*t1y-t3y)/(t3y-t0y);
ssy2=(2*t2y-t3y)/(t3y-t0y);
T1y=[-t3y/2:dty:+t3y/2];
daty1=length(T1y);
Ty=[0:dty:t4y];
daty=length(Ty);
Vmaxyf=120;
Vmaxyb=74;
Volyf=zeros(1,daty);  %Forward
Volyb=zeros(1,daty);  %Backward
Volyf(1:daty1)=Vmaxyf*tripuls(T1y,t3y,ssy1);  %Forward
Volyb(1:daty1)=Vmaxyb*tripuls(T1y,t3y,ssy2);  %Backward
m1x=400e-3;
m2x=300e-3; 
m1y=(m1x+m2x+Mp)+96e-3;
m2y=200e-3; 
MUky=0.51;
MUsy=0.61;
Fny=(m1y+m2y+Mp)*g;
Fsy=MUsy*Fny;                      
Fcy=MUky*Fny;
%% ------------------------------------- Controller Gains
% Controller Gains
Kpro=1/sqrt(0.73)*1e3;
Kint=0;
Kder=0;
%% ------------------------------------- Initial Condition
% X
hx=zeros(1,datx);
Fhx=zeros(1,datx);
sx=ones(1,datx)*eps*1e6;
x1(1,1)=0;
x1(1,2)=eps;
x2(1)=0;
x2(2)=eps;
% Y
hy=zeros(1,daty);
Fhy=zeros(1,daty);
sy=ones(1,daty)*eps*1e6;
y1(1,1)=0;
y1(1,2)=eps;
y2(1)=0;
y2(2)=eps;

%% ------------------------------------- System Equeitions
% ------------------------- X Started
tic
kx=1
xf(1)=0;
%Control
errorx(kx)=Xd-xf(1);
Uc=Kpro*sqrt(abs(errorx(kx)))*sign(errorx(kx));%Kpro*sqrt(errorx(kx));
if abs(Uc)>1
    if errorx(kx)>=0
    Volx(1,:)=Volxf(1,:);    
    else
    Volx(1,:)=Volxb(1,:);    
    end
else
    if errorx(kx)>=0
    Volx(1,:)=abs(Uc)*Volxf(1,:); 
    else
    Volx(1,:)=abs(Uc)*Volxb(1,:);
    end
end
MAL=0.1e-6;
while abs(errorx(kx))>MAL
% Puls
for i=2:1:datx
dVolx(i)=(Volx(kx,i)-Volx(kx,i-1))/dtx;
% piezo force
hx(i)=(ALFp*de*dVolx(i)-BETp*abs(dVolx(i))*hx(i-1)-GAMp*dVolx(i)*abs(hx(i-1)))*dtx+hx(i-1);
Fpx(i)=kp*(de*Volx(kx,i)-hx(i));
% friction
dx1(i)=(x1(kx,i)-x1(kx,i-1))/dtx;
sx(i)=sign(dx1(i))*(Fcx+(Fsx-Fcx)*exp(-(abs(dx1(i))/vs)^DEL));     % constant velocity behavior
dzx(i)=dx1(i)*(1-sign(Fhx(i)/sx(i))*(abs(Fhx(i)/sx(i)))^n);
Fhx(i+1)=(ALF*dzx(i)-BET*abs(dzx(i))*Fhx(i)-GAM*dzx(i)*abs(Fhx(i)))*dtx+Fhx(i);       % hysteresis friction
Ffx(i)=Fhx(i)+SIG1*dzx(i)+SIG2*dx1(i);                           % Frictional force
Ffx(i)=abs(Ffx(i))*sign(dx1(i));
% governing Equation
x1(kx,i+1)=2*x1(kx,i)-x1(kx,i-1)+(-kp*(x1(kx,i)-x2(kx,i))-cp*(x1(kx,i)-x1(kx,i-1)-x2(kx,i)+x2(kx,i-1))/dtx+Fpx(i)-Ffx(i))*dtx^2/m1x;
x2(kx,i+1)=2*x2(kx,i)-x2(kx,i-1)+(+kp*(x1(kx,i)-x2(kx,i))+cp*(x1(kx,i)-x1(kx,i-1)-x2(kx,i)+x2(kx,i-1))/dtx-Fpx(i))*dtx^2/m2x;
end
TTx(kx,:)=Tx+(kx-1)*t4x;
% Control
xf(kx+1)=x1(kx,datx)+xf(kx);
dxf(kx+1)=(xf(kx+1)-xf(kx))/t4x;
errorx(kx+1)=Xd-xf(kx+1);
Uc=Kpro*sqrt(abs(errorx(kx+1)))*sign(errorx(kx+1));%sqrt(errorx(kx+1));
if abs(Uc)>1
    if errorx(kx+1)>=0
    Volx(kx+1,:)=Volxf(1,:);    
    else
    Volx(kx+1,:)=Volxb(1,:);    
    end
else
    if errorx(kx+1)>=0
    Volx(kx+1,:)=abs(Uc)*Volxf(1,:); 
    else
    Volx(kx+1,:)=abs(Uc)*Volxb(1,:);
    end
end
% Next Step
hx=zeros(1,datx);
Fhx=zeros(1,datx);
sx=ones(1,datx)*eps*1e6;
x1(kx+1,1)=0;
x1(kx+1,2)=eps;
x2(kx+1,1)=0;
x2(kx+1,2)=eps;
% next kx
kx=kx+1
end
% ------------------------- X Finished
% ------------------------- Y Started
% y
ky=1
yf(1)=0;
%Control
errory(ky)=Yd-yf(1);
Uc=Kpro*sqrt(abs(errorx(ky)))*sign(errory(ky));%sqrt(errory(ky));
if abs(Uc)>1
    if errory(ky)>=0
    Voly(1,:)=Volyf(1,:);    
    else
    Voly(1,:)=Volyb(1,:);    
    end
else
    if errory(ky)>=0
    Voly(1,:)=abs(Uc)*Volyf(1,:); 
    else
    Voly(1,:)=abs(Uc)*Volyb(1,:);
    end
end
while abs(errory(ky))>MAL 
% Puls
for i=2:1:daty
dVoly(i)=(Voly(ky,i)-Voly(ky,i-1))/dty;
% piezo force
hy(i)=(ALFp*de*dVoly(i)-BETp*abs(dVoly(i))*hy(i-1)-GAMp*dVoly(i)*abs(hy(i-1)))*dty+hy(i-1);
Fpy(i)=kp*(de*Voly(ky,i)-hy(i));
% friction
dy1(i)=(y1(ky,i)-y1(ky,i-1))/dty;
sy(i)=sign(dy1(i))*(Fcy+(Fsy-Fcy)*exp(-(abs(dy1(i))/vs)^DEL));     % constant velocity behavior
dzy(i)=dy1(i)*(1-sign(Fhy(i)/sy(i))*(abs(Fhy(i)/sy(i)))^n);
Fhy(i+1)=(ALF*dzy(i)-BET*abs(dzy(i))*Fhy(i)-GAM*dzy(i)*abs(Fhy(i)))*dty+Fhy(i);       % hysteresis friction
Ffy(i)=Fhy(i)+SIG1*dzy(i)+SIG2*dy1(i);                           % Frictional force
Ffy(i)=abs(Ffy(i))*sign(dy1(i));
% governing Equation
y1(ky,i+1)=2*y1(ky,i)-y1(ky,i-1)+(-kp*(y1(ky,i)-y2(ky,i))-cp*(y1(ky,i)-y1(ky,i-1)-y2(ky,i)+y2(ky,i-1))/dty+Fpy(i)-Ffy(i))*dty^2/m1y;
y2(ky,i+1)=2*y2(ky,i)-y2(ky,i-1)+(+kp*(y1(ky,i)-y2(ky,i))+cp*(y1(ky,i)-y1(ky,i-1)-y2(ky,i)+y2(ky,i-1))/dty-Fpy(i))*dty^2/m2y;
end
TTy(ky,:)=Ty+(ky-1)*t4y;
% Control
yf(ky+1)=y1(ky,daty)+yf(ky);
dyf(ky+1)=(yf(ky+1)-yf(ky))/t4y;
errory(ky+1)=Yd-yf(ky+1);
Uc=Kpro*sqrt(abs(errory(ky+1)))*sign(errory(ky+1));%sqrt(errory(ky+1));
if abs(Uc)>1
    if errory(ky+1)>=0
    Voly(ky+1,:)=Volyf(1,:);    
    else
    Voly(ky+1,:)=Volyb(1,:);    
    end
else
    if errory(ky+1)>=0
    Voly(ky+1,:)=abs(Uc)*Volyf(1,:); 
    else
    Voly(ky+1,:)=abs(Uc)*Volyb(1,:);
    end
end
% Next Step
hy=zeros(1,daty);
Fhy=zeros(1,daty);
sy=ones(1,daty)*eps*1e6;
y1(ky+1,1)=0;
y1(ky+1,2)=eps;
y2(ky+1,1)=0;
y2(ky+1,2)=eps;
% next ky
ky=ky+1
end
% ------------------------- Y Finished
toc
%% ------------------------------------- Classification
% X
x1(kx,:)=[];
x2(kx,:)=[];
j=0;
for j=1:1:kx-1
    xx1(1,datx*(j-1)+1:datx*j)=x1(j,1:datx)+xf(j);
    TTx1(1,datx*(j-1)+1:datx*j)=TTx(j,1:datx);
    Volxx(1,datx*(j-1)+1:datx*j)=Volx(j,1:datx);
end
% Y
y1(ky,:)=[];
y2(ky,:)=[];
j=0;
for j=1:1:ky-1
    yy1(1,daty*(j-1)+1:daty*j)=y1(j,1:daty)+yf(j);
    TTy1(1,daty*(j-1)+1:daty*j)=TTy(j,1:daty);
    Volyy(1,daty*(j-1)+1:daty*j)=Voly(j,1:daty);
end

%% ------------------------------------- Result
AMXP=(Xd+MAL)*ones(1,length(TTx1));
AMXN=(Xd-MAL)*ones(1,length(TTx1));
AMYP=(Yd+MAL)*ones(1,length(TTy1));
AMYN=(Yd-MAL)*ones(1,length(TTy1));


% X
figure(1)
subplot(2,2,1)
plot(TTx1*1e3,Volxx*max(abs(xx1*1e6))/Vmaxxb,'k-.','linewidth',1.8)
hold on
plot(TTx1*1e3,xx1*1e6,'linewidth',2.5,'color',[1,0,0])
hold on
plot(TTx1*1e3,Xd*ones(1,length(TTx1))*1e6,'b--','linewidth',1.5)
hold on
plot(TTx1*1e3,AMXP*1e6,'g-','linewidth',1.5)
hold on
plot(TTx1*1e3,AMXN*1e6,'g-','linewidth',1.5)
xlabel('Time (ms)','fontsize',14)
ylabel('Displacement (\mum)','fontsize',14)
title('X-Square Root','fontsize',18)
%grid on
box on


% Y
%figure(2)
subplot(2,2,3)
plot(TTy1*1e3,Volyy*max(abs(yy1*1e6))/Vmaxyb,'k-.','linewidth',1.8)
hold on
plot(TTy1*1e3,yy1*1e6,'linewidth',2.5,'color',[1,0,0])
hold on
plot(TTy1*1e3,Yd*ones(1,length(TTy1))*1e6,'b--','linewidth',1.5)
hold on
plot(TTy1*1e3,AMYP*1e6,'g-','linewidth',1.5)
hold on
plot(TTy1*1e3,AMYN*1e6,'g-','linewidth',1.5)
xlabel('Time (ms)','fontsize',14)
ylabel('Displacement (\mum)','fontsize',14)
title('Y-Square Root','fontsize',18)
%grid on
box on