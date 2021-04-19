clc;clear%%�����������������ʼ��������

%% Model parametersģ�ͳ�ʼ������������ѧϵͳ�����ĳ�ʼֵ��
c=0;beta=0;

delta_I0=0.13; delta_q=0.13;

gama_I=0.007; gama_H=0.014;

q=0; alpha=2.7e-4;

theta=1; lam=1/14; k=0;%�������߸�����

T=30; t=0.1; NN=T/t;%%����Ԥ������Ϊ60�죬��0.02Ϊ�������ODE



%% Initial values����ѧ��������ʼֵ
%% Initial values����ѧ��������ʼֵ
S=20000000; E=27; I=161; Sq=0;

Eq=0; H=I+Eq; R=1880; D=7; sigma=1/5.5;

AA=[S E I Sq Eq H R D];

for ii=1:NN

%% Increased isolation speed of patient due to new
%% hospitals���ܽ����人С��ɽҽԺʹ�ø����ٶ�������ǰ�����ܺ���в��

if (ii*t)>=14

delta_I=delta_I0*0;

else

delta_I=delta_I0;

end

%% Modified SEIR Transmission dynamics model  �Ľ���SEIR��Ⱦ��ģ�Ͷ���ѧ����

dS = -(beta*c+c*q*(1-beta))*S*(I+theta*E)+lam*Sq;

dE = beta*c*(1-q)*S*(I+theta*E)-sigma*E;

dI = sigma*E-(delta_I+alpha+gama_I)*I+k*R;

dSq = (1-beta)*c*q*S*(I+theta*E)-lam*Sq;

dEq = beta*c*q*S*(I+theta*E)-delta_q*Eq;

dH = delta_I*I+delta_q*Eq-(alpha+gama_H)*H;

dR = gama_I*I+gama_H*H-k*R;

dD = alpha*(I+H);

%% Euler integration algorithm  ODE��ŷ����ֵ��

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
%Syi(:,1)= round(AA(1:1/t:size(AA,1),1));%��ǰ�׸�������

%Eqian(:,1)= round(AA(1:1/t:size(AA,1),2));%��ǰǱ��������

Infected(:,1)= round(AA(1:1/t:size(AA,1),3));%��ǰ��Ⱦ������===============

Cured(:,1)=round(AA(1:1/t:size(AA,1),7));%��ǰ����������====================

Death(:,1)= round(AA(1:1/t:size(AA,1),8));%��ǰ��������=================

All(:,1)=Death(:,1);%Infected(:,1)+Cured(:,1)+Death(:,1);


%%plot graph  ����ͼ��
plot(0:T,[Death ]);grid on;%Syi Cured Death
xlabel('��');ylabel('����')
legend( '��ǰ��������')%'��ǰ�׸�������',,'��ǰ��������', '��ǰ��������'
title('ģ������Ԥ�������')
