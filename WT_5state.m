clear; close all;clc;
%% 1.5MW改8MW基本参数
N_B=3;%叶片数
Ng=97;%齿轮箱转速比
rou=1.225;%空气密度
R=63;%转子半径
Jr=35444067;%风轮转动惯量
Jg=534.116;%发电机惯量
Ks=867637000;%传动系扭转刚度
Ds=6215000;%传动扭转阻尼系数
tau_beta =0.05;%桨距动作时间常数
tau_g=0.1;%发电机时间常数
eta_g=0.944;%发电机效率
Wr_rated =12.1*pi/30;%转子稳态转速
Wg_rated=1173.7*pi/30;%发电机转速
Tg_rated =43093.55;%转矩
h=87.6;%塔架高度
theta_0=Tg_rated*Ng/Ks;
Wr_0=Wr_rated;
Wg_0=Wg_rated;
%% 转矩，推力偏导数
syms Wr beta V
%CP Ct非线性公式
lamda = Wr*R/V;%叶尖速比
lamda_i=1/(lamda+0.08*beta)-0.035/(beta^3+1);%中间变量
Cp=6.909*(7.022*lamda_i-0.04176*beta-0.3863)*exp(-14.52*lamda_i); % 功率系数
% Ct=Cp;%这条公式有问题需要去找Ct关系式 为了先写C矩阵
% Cp=0.5176*(116*lamda_i-0.4*beta-5)*exp(-21*lamda_i)+0.0068*lamda;
Ct=0.08698-0.003371*beta-0.05895*lamda+0.00567*beta*lamda+0.06499*lamda^2-0.003096*beta*lamda^2-0.009577*lamda^3+0.0002103*beta*lamda^3+0.0005667*lamda^4-(5.24*10^-6)*beta*lamda^4-(1.199*10^-5)*lamda^5;
Pr(Wr,beta,V)=0.5*rou*pi*R^2*Cp*V^3;% wind power captured by the rotor [w]
Tr(Wr,beta,V)=Pr(Wr,beta,V)/Wr;%aerodynamic torque applied to the rotor [Nm]
Ft(Wr,beta,V)=0.5*rou*pi*R^2*V^2*Ct;%thrust force applied to the rotor [N]
%计算稳态工况点
Tr_rated=43093.55*Ng;%额定转子转矩
Wr_Rated=12.1*pi/30;% 额定转子转速从转每分到弧度每秒
Pr_rated=5e6;%额定功率
V_op=1:1:25;
beta_op=zeros(0,25);% Initialize the pitch to zero
lamda_op =7.55*ones(1,11);
for Vi=12:1:25
    sol =vpasolve(Pr == Pr_rated, Wr==Wr_Rated, V==Vi, [Wr,beta,V]);
%     disp(so1.beta)
%     disp(so1.V)
    V_op(Vi)=sol.V;
    beta_op(Vi)=sol.beta;
    lamda_op(Vi)=Wr_rated*R/sol.V;
end
%画稳态工况轨迹
figure(1)
subplot(2,1,1)
plot(V_op,beta_op,'LineWidth',1.2)
xlim([8,25]);
%xlabel("Wind Speed(m/s)","Fontname",Arial')
%ylabel('β(deg)','Fontname',"Arial')
xlabel('风速')
ylabel('变桨角')
subplot(2,1,2)
plot(V_op,lamda_op,'LineWidth',1.2)
xlim([8,25]);
xlabel('风速')
ylabel('叶尖速比')
%稳态工况点处Tr的偏导曲线
K_Tw=diff(Tr,'Wr',1);
K_Tbeta=diff(Tr,'beta',1);
K_TV=diff(Tr,'V',1);
KK_Tw =zeros(1,25);
KK_Tbeta=zeros(1,25);
KK_TV=zeros(1,25);
for i=12:1:25
    K_Tw_value = double(subs(K_Tw, {Wr beta V}, {Wr_rated, beta_op(i),V_op(i)}));
    K_Tbeta_value = double(subs(K_Tbeta, {Wr beta V}, {Wr_rated, beta_op(i),V_op(i)}));
    K_TV_value =double(subs(K_TV, {Wr beta V},{Wr_rated, beta_op(i),V_op(i)}));
    KK_Tw(i)=K_Tw_value;
    KK_Tbeta(i)=K_Tbeta_value;
    KK_TV(i)= K_TV_value;
end
K_Fw=diff(Ft,'Wr',1);
K_Fbeta=diff(Ft,'beta',1);
K_FV=diff(Ft,'V',1);
KK_Fw =zeros(1,25);
KK_Fbeta=zeros(1,25);
KK_FV=zeros(1,25);
for i=12:1:25
    K_Fw_value = double(subs(K_Fw, {Wr beta V}, {Wr_rated, beta_op(i),V_op(i)}));
    K_Fbeta_value = double(subs(K_Fbeta, {Wr beta V}, {Wr_rated, beta_op(i),V_op(i)}));
    K_FV_value =double(subs(K_FV, {Wr beta V},{Wr_rated, beta_op(i),V_op(i)}));
    KK_Fw(i)=K_Fw_value;
    KK_Fbeta(i)=K_Fbeta_value;
    KK_FV(i)= K_FV_value;
end
figure(2)
subplot(3,1,1)
plot(V_op,KK_Tw,'LineWidth',1.2)
xlim([8,25]);
xlabel('风速')
ylabel('Tr对Wr求导')
subplot(3,1,2)
plot(V_op,KK_Tbeta,'LineWidth',1.2)
xlim([8,25]);
xlabel('风速')
ylabel('Tr对beta求导')
subplot(3,1,3)
plot(V_op,KK_TV,'LineWidth',1.2)
xlim([8,25]);
xlabel('风速')
ylabel('Tr对v求导')
for Vi=12:1:25
    A=[ 0            1               -1/Ng            0              0   ;
      -Ks/Jr     (KK_Tw(Vi)-Ds)/Jr   Ds/(Jr*Ng)   KK_Tbeta(Vi)/Jr    0   ;
      Ks/(Jg*Ng)     Ds/(Jg*Ng)      -Ds/(Jg*Ng^2)    0            -1/Jg ;
        0            0                0           -1/tau_beta        0   ;
        0            0                0                 0       -1/tau_g ;
      ];
    B=[0 0;0 0;0 0;1/tau_beta 0;0 1/tau_g]; 
    Bv=[0;KK_TV(Vi)/Jr;0;0;0];
%     C=[0,0,1,0,0];
    C=[0  0  eta_g*Tg_rated 0 eta_g*Wg_rated
      h*KK_Fw(Vi) 0 0 0 h*KK_Fbeta(Vi)
      Ks Ds -Ds/Ng  0  0];
%        Ctomigar 0 0 0 Ctbeta];
%     D=zeros(3,2);
    D=[0];
    Wv=[0 h*KK_FV(i) 0];
states ={'theta','Wr','We','beta','Te'};
inputs ={'beta_ref','Te_ref','V'};
outputs ={'Pe Mt Ts Ct'};
sys =ss(A,[B Bv],C,D,'statename',states,"inputname",inputs,'outputname',outputs);
plant(:,:,Vi)= ss(A,[B Bv],C,D,'statename',states,'inputname',inputs,'outputname',outputs);
plant(:,:,Vi,1)
end
A_at_12m_s = plant(:,:,12,1);
disp('A矩阵（风速为12 m/s）：');
disp(A_at_12m_s);
% eig(A);
% Qc=ctrb(A,B);
% if rank(Qc)==rank(A)
%     disp('System is controllable');
% else
%     disp('System is not controllable');
% end
% sys=ss(A,B,C,D)