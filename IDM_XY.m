% In the name of Allah
% Ali Abdi
% 1397-09-22
% Thesis
clear all
close all
clc
% Main Code
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
Vmaxxf=100
%Vmaxxb=49;
Volxf=zeros(1,datx);  %Forward
Volxb=zeros(1,datx);  %Backward
Volxf(1:datx1)=Vmaxxf*tripuls(T1x,t3x,ssx1);  %Forward
%Volxb(1:datx1)=Vmaxxb*tripuls(T1x,t3x,ssx2);  %Backward
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
kx=1
xf(1)=0;
Volx(1,:)=Volxf(1,:); 

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

%% ------------------------------------- Result
% X
figure(1)
%subplot(2,2,1)
plot(TTx*1e3,Volx*max(abs(x1*1e6))/Vmaxxf,'k-.','linewidth',1.8)
hold on
plot(TTx*1e3,x1(1:length(TTx))*1e6,'linewidth',2.5,'color',[1,0,0])

xlabel('Time (ms)','fontsize',14)
ylabel('Displacement (\mum)','fontsize',14)
title('X-Cascade','fontsize',18)
grid on
