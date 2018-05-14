%%%%һά����ά�����㷨
clc;
clear all;
close all;
%% ��֪����
c=3e8;%����
R_ref=0;%��ο�Ŀ��λ��
x_target=0;
y_target=0;%ʵ��Ŀ�������λ��
fc=183e9;%�����źų�ʼƵ��
B=150e8;%�����
Tr=1.5e-6;%������
gamma=B/Tr;%��Ƶб��
%% �趨�շ����ߵ�λ��
x_TR=0;
y_TR=0;
z_TR=0;
%% ȷ������Ƶ���λ��
Nf=256;%ɨ���Ƶ����
% f_step=B/Nf;%Ƶ�ʱ仯���
% f=(0:Nf-1)*f_step;%���Ե�Ƶ�ź�
% f=linspace(0,B,Nf);
% z_d=c*t/2;
%% ����ά����ķֱ��� 
t=linspace(0,Tr,Nf);
z_d=linspace(0,Tr*c/2,Nf);
f=gamma*t;
delta_z=c/(2*B);
z_grid=[(0:Nf-1)]*delta_z;
%% �趨Ŀ����λ��
Ntar=8;%Ŀ��ĸ���
Ptar=[1*delta_z 1
    20*delta_z 1
    30*delta_z 1
    40*delta_z 1
    50*delta_z 1
    60*delta_z 1
    70*delta_z 1
    200*delta_z 1];
%% �ز�����
s_if=0;%��ʼ����Ƶ�ź�  
for i=1:Ntar;
z_target=Ptar(:,1);
sigma=Ptar(:,2);
R_i=sqrt((x_TR-x_target)^2+(y_TR-y_target)^2+(z_TR-z_target(i)).^2);%Ŀ�굽�շ�����֮��ľ���
phase_r=2*pi*(fc*(t-(2*R_i/c))+gamma*((t-(2*R_i/c)).^2)/2);
s_r=sigma(i)*exp(-1j*phase_r);
phase_ref=2*pi*(fc*(t-(2*R_ref/c))+gamma*((t-(2*R_ref/c)).^2)/2);
s_ref=exp(-1j*phase_ref);
s_if=s_if+s_r.*conj(s_ref);%���в�Ƶ����
end         
% figure
% plot(fftshift(abs(fft(s_if))));
%% �Ƚ���FFT����
s_if=fft(s_if);
%% RVP���
s_comp=s_if.*exp(-1j*pi*f.^2/gamma);
% figure
% plot(abs((ifft((s_comp)))));   
%% ���ο��źŵ���λ�������
k_i=2*pi*(f+fc)/c;
s_dref=s_comp.*exp(1j*2*k_i*R_ref);
%% ȷ��������Ƶ��
kx=abs(k_i*(x_TR-x_target)/R_i);
ky=abs(k_i*(y_TR-y_target)/R_i);
kz=(sqrt(4*(k_i.^2)-kx.^2-ky.^2));
% ky=0;
% kx=0;
% kz=2*k_i;
%% �Է�λά���в���
s_compa=s_dref.*exp(-1j*(kx*(x_TR-x_target)+ky*(y_TR-y_target)));
%% �Ծ���ά���в���
% s_compd=s_compa.*exp(-1j*kz*z_TR);
%% �Ծ���ά���и���Ҷ��任
% G=s_compd.*exp(1j*kz.*z_d);
G=ifft(s_compa);
%% ��ѹ�ӹ�һ������
G_pc=fft(G)./max(fft(G));
figure
plot(z_grid,(abs(G_pc)),Ptar(:,1),Ptar(:,2),'*');
title('һά����ά�������ͼ');
xlabel('����/m');
ylabel('��һ������');
