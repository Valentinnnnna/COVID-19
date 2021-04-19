clc;clear%%清理环境变量，构造初始运行条件

%% Model parameters模型初始参数（即动力学系统参数的初始值）
c=0;beta=0;

delta_I0=0.13; delta_q=0.13;

gama_I=0.007; gama_H=0.014;

q=0; alpha=2.7e-4;

theta=1; lam=1/14; k=0;%治愈患者复阳率

T=30; t=0.1; NN=T/t;%%设置预测天数为60天，以0.02为步长求解ODE



%% Initial values动力学各变量初始值
%% Initial values动力学各变量初始值
S=20000000; E=27; I=161; Sq=0;

Eq=0; H=I+Eq; R=1880; D=7; sigma=1/5.5;

AA=[S E I Sq Eq H R D];

for ii=1:NN

%% Increased isolation speed of patient due to new
%% hospitals两周建成武汉小汤山医院使得隔离速度在两周前和两周后具有差别

if (ii*t)>=14

delta_I=delta_I0*0;

else

delta_I=delta_I0;

end

%% Modified SEIR Transmission dynamics model  改进的SEIR传染病模型动力学方程

dS = -(beta*c+c*q*(1-beta))*S*(I+theta*E)+lam*Sq;

dE = beta*c*(1-q)*S*(I+theta*E)-sigma*E;

dI = sigma*E-(delta_I+alpha+gama_I)*I+k*R;

dSq = (1-beta)*c*q*S*(I+theta*E)-lam*Sq;

dEq = beta*c*q*S*(I+theta*E)-delta_q*Eq;

dH = delta_I*I+delta_q*Eq-(alpha+gama_H)*H;

dR = gama_I*I+gama_H*H-k*R;

dD = alpha*(I+H);

%% Euler integration algorithm  ODE的欧拉数值解

S =S+dS*t;

E = E+dE*t;

I = I+dI*t;

Sq = Sq+dSq*t;

Eq = Eq+dEq*t;

H = H+dH*t;

R = R+dR*t;

D = D+dD*t;

AA=[AA; S E I Sq Eq H R D];

end

%% Theoretical estimation
%AA=[S E I Sq Eq H R D];  Sq Eq H
%Syi(:,1)= round(AA(1:1/t:size(AA,1),1));%当前易感者人数

%Eqian(:,1)= round(AA(1:1/t:size(AA,1),2));%当前潜伏者人数

Infected(:,1)= round(AA(1:1/t:size(AA,1),3));%当前感染者人数===============

Cured(:,1)=round(AA(1:1/t:size(AA,1),7));%当前治愈者人数====================

Death(:,1)= round(AA(1:1/t:size(AA,1),8));%当前死亡人数=================

All(:,1)=Death(:,1);%Infected(:,1)+Cured(:,1)+Death(:,1);


%%plot graph  绘制图像
plot(0:T,[Death ]);grid on;%Syi Cured Death
xlabel('天');ylabel('人数')
legend( '当前死亡人数')%'当前易感者人数',,'当前治愈人数', '当前死亡人数'
title('模型疫情预测结果情况')
