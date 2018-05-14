%%%%一维距离维成像算法
clc;
clear all;
close all;
%% 已知条件
c=3e8;%光速
R_ref=0;%设参考目标位置
x_target=0;
y_target=0;%实际目标的坐标位置
fc=183e9;%发射信号初始频率
B=150e8;%设带宽
Tr=1.5e-6;%设脉宽
gamma=B/Tr;%调频斜率
%% 设定收发天线的位置
x_TR=0;
y_TR=0;
z_TR=0;
%% 确定采样频点的位置
Nf=256;%扫描的频点数
% f_step=B/Nf;%频率变化间隔
% f=(0:Nf-1)*f_step;%线性调频信号
% f=linspace(0,B,Nf);
% z_d=c*t/2;
%% 距离维方向的分辨率 
t=linspace(0,Tr,Nf);
z_d=linspace(0,Tr*c/2,Nf);
f=gamma*t;
delta_z=c/(2*B);
z_grid=[(0:Nf-1)]*delta_z;
%% 设定目标点的位置
Ntar=8;%目标的个数
Ptar=[1*delta_z 1
    20*delta_z 1
    30*delta_z 1
    40*delta_z 1
    50*delta_z 1
    60*delta_z 1
    70*delta_z 1
    200*delta_z 1];
%% 回波构建
s_if=0;%初始化中频信号  
for i=1:Ntar;
z_target=Ptar(:,1);
sigma=Ptar(:,2);
R_i=sqrt((x_TR-x_target)^2+(y_TR-y_target)^2+(z_TR-z_target(i)).^2);%目标到收发天线之间的距离
phase_r=2*pi*(fc*(t-(2*R_i/c))+gamma*((t-(2*R_i/c)).^2)/2);
s_r=sigma(i)*exp(-1j*phase_r);
phase_ref=2*pi*(fc*(t-(2*R_ref/c))+gamma*((t-(2*R_ref/c)).^2)/2);
s_ref=exp(-1j*phase_ref);
s_if=s_if+s_r.*conj(s_ref);%进行差频处理
end         
% figure
% plot(fftshift(abs(fft(s_if))));
%% 先进行FFT处理
s_if=fft(s_if);
%% RVP项补偿
s_comp=s_if.*exp(-1j*pi*f.^2/gamma);
% figure
% plot(abs((ifft((s_comp)))));   
%% 将参考信号的相位误差消除
k_i=2*pi*(f+fc)/c;
s_dref=s_comp.*exp(1j*2*k_i*R_ref);
%% 确定波数域频率
kx=abs(k_i*(x_TR-x_target)/R_i);
ky=abs(k_i*(y_TR-y_target)/R_i);
kz=(sqrt(4*(k_i.^2)-kx.^2-ky.^2));
% ky=0;
% kx=0;
% kz=2*k_i;
%% 对方位维进行补偿
s_compa=s_dref.*exp(-1j*(kx*(x_TR-x_target)+ky*(y_TR-y_target)));
%% 对距离维进行补偿
% s_compd=s_compa.*exp(-1j*kz*z_TR);
%% 对距离维进行傅里叶逆变换
% G=s_compd.*exp(1j*kz.*z_d);
G=ifft(s_compa);
%% 脉压加归一化处理
G_pc=fft(G)./max(fft(G));
figure
plot(z_grid,(abs(G_pc)),Ptar(:,1),Ptar(:,2),'*');
title('一维距离维成像仿真图');
xlabel('距离/m');
ylabel('归一化幅度');
