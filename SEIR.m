clc;clear%%�����������������ʼ��������

%% Model parametersģ�ͳ�ʼ������������ѧϵͳ�����ĳ�ʼֵ��
c=3;beta=10.9e-9;%cΪ�׸���Ⱥƽ��ÿ�սӴ�����,betaΪ��Ⱦ����

delta_I0=0.13; delta_q=0.13;%delta_I0Ϊ��Ⱦ��Ⱥ�������������ʣ�delta_qΪ����Ǳ����Ⱥ��������������

gama_I=0.007; gama_H=0.014;% gama_IΪ��Ⱦ��Ⱥ�ָ���,gama_HΪ����������Ⱥ�Ļָ�����

q=1e-6; alpha=2.7e-4;% qΪ�׸���Ⱥ����۲����,alphaΪ������

theta=1; lam=1/14; k=0;%thetaΪǱ����Ⱥ��Ⱦ�׸���Ⱥ���ʣ�lamdaΪ����۲���Ⱥ�����۲�����,k�������߸�����


T=30; t=0.1; NN=T/t;%%����Ԥ������Ϊ30�죬��0.1Ϊ�������ODE



%% Initial values����ѧ��������ʼֵ
%% Initial values����ѧ��������ʼֵ
S=3.282e8; E=288260; I=3973973; Sq=(28778117-28030005)*c;

Eq=Sq*400/2776; H=I+Eq; R=24346766; D=529759; sigma=1/5.5;

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

%Cured(:,1)=round(AA(1:1/t:size(AA,1),7));%��ǰ����������====================

Death(:,1)= round(AA(1:1/t:size(AA,1),8));%��ǰ��������=================

All(:,1)=Infected(:,1)%+Cured(:,1)+Death(:,1);


%%plot graph  ����ͼ��
plot(0:T,[Death ]);grid on;%Syi Cured Death
xlabel('��');ylabel('����')
legend( '��ǰ��������')%'��ǰ�׸�������','��ǰ��Ⱦ����','��ǰ��������', '��ǰ��������'
title('ģ������Ԥ�������')
%% raw data  ��ʵ�е���ʵȷ���������ں�ģ��Ԥ��Ľ�����Աȣ�����ģ��Ч����˵���������Ժ��ģ��Ԥ����������ʵ��

hold on

z_Infected=[3973973 3958455 3945678 3922861 3922679 3816399 3705295 3710142 3598717]';%�ִ�ȷ��

z_Cured=[24346766 24423589 24480522 24560856 24626410 24796161 24970980 25236187 25305332]';%�ۼ�����

z_Dead=[529759 530526 531818 534211 537076 542141 543740 547589 555293]';%�ۼ�����


plot([1:length(z_Infected)]'-1,[z_Dead ],'*')%z_Cured z_Dead