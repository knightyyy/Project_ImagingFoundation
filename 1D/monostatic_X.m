%% 一维方向维成像算法
clc;
clear all;
close all;
%% 已知条件
c=3e8;%光速
fc=183e9;%雷达发射信号中心频率
Yc=0.5;%目标与探测器移动平面的距离0.5m
v=0.001;%探测器移动速度1mm/s
L=0.1;%合成孔径一半的距离10cm
theta=10*pi/180;%天线方向角
%% 合成孔径的参数设定
lambda=c/fc;%波长
w=2*pi*fc;
k=2*w/c;%波数域频率
D=2*L;%孔径实际大小
%% 设定目标点所在位置的有效范围
X_range=2*L+2*Yc*tan(theta/2);
%% 方位向分辨率
delta_x=Yc*lambda/(2*D);
% Na=2*L/delta_x;%采样点的个数
Na=256;
x=linspace(-L,L,Na);%探测器运动范围
u=x/v;
% du=2*L/v/Na;
% fu=linspace(-1/2/du,1/2/du,Na);%fu域序列
%% 目标区域范围
xa=linspace(-X_range/2,X_range/2,Na);
%% 设定方位向距离及反射系数
% Ntar=2;
% Ptar=[0,1
%     4*delta_x,1
%     10*delta_x,1
%     20*delta_x,1
%     40*delta_x,1
%     ];
Ntar=2;
G=X_range/2/delta_x;
Ptar=[G*rand(1,Ntar)*delta_x].';%随机设置目标点的位置
% Ntar=20;
% Ptar=linspace(-10*delta_x,10*delta_x,Ntar);
%% 构建回波信号
s=0;
for i=1:1:Ntar; 
%     xn=Ptar(i,1);
%     sigma=Ptar(i,2);
    xn=Ptar(i);
    sigma=1;
    R=sqrt((Yc^2)+((xn-x).^2));
    s=s+sigma*exp(-1j*k*R);%回波信号 如果目标移动范围大于L,小于X_range/2时，在-L,L的范围内无法表示出来
end
% figure
% plot(abs(fft(s)));
% SNR_db=10; 
% SNR=10^(SNR_db/10); 
% %%%%%%% %计算信号能量 
% xt=s;
% sigpower=0; 
% for i=1:length(xt) 
%     sigpower=sigpower+xt(i).^2;
% end
% %生成噪声 
% noisepower=sigpower/SNR;%所需噪声能量 
% noise=randn(1,length(xt));%产生均值为零，方差为一的随机高斯序列 
% noise=noise-repmat(mean(mean(noise)),1,length(xt));%让均值更接近于零%噪声校正 
% AWGNpower=0; 
% for i=1:length(xt) 
%     AWGNpower=AWGNpower+noise(1,length(xt)).^2;%生成噪声的能量
% end
% K=noisepower/AWGNpower;%所需噪声能量与实际生成噪声能量的比值 
% K1=sqrt(K);%校正生成噪声加权值 
% snr=10*log10(sigpower/sum((K1*noise).^2));%实际输出信噪比
% s=xt+(K1*noise);
%% 对回波在慢时间域内做傅里叶变换
% Kx=(k*((x-xa))./sqrt((Yc^2)+(((x-xa)).^2)));
kx_range=2*pi*X_range/(lambda*sqrt((X_range/2).^2+Yc.^2));
Kx=linspace(-kx_range,kx_range,Na);
Sn=exp(-1j*Kx.'*x)*s.';
Sn=Sn.';
% Sn=fftshift(fft(fftshift(s)));
%% 进行相位补偿
%     Kx=(k*(xn-x)./sqrt((Yc^2)+((xn-x).^2)));
val=exp(1j*((sqrt(((k.^2))-(Kx.^2))))*Yc);                                                                                                                                                                                                                                                                                                                                                                                            
S_comp=Sn.*val;
%% 对方向位进行反傅里叶变换                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
% s_ifft=fftshift(ifft(fftshift(S_comp)));
s_ifft=exp(1j*xa.'*Kx)*S_comp.';
s_ifft=s_ifft.';
%% 对幅值进行归一化处理
s_ifft=s_ifft./max(s_ifft);
% %% 画图
% figure                                      
% plot(xa,abs(s_ifft));
% xlabel('位置/m');
% ylabel('幅值');
% axis on;
% s_fft=fftshift(fft(fftshift(s_ifft)));
figure
plot(xa,abs(s_ifft));
hold on
plot(Ptar(:),ones(1,length(Ptar)),'r*');
% title('频域图');
% vm=sqrt(4-((D.^2)/((L.^2)+Yc.^2)));
% delta_z=c/(2*(2-vm)*fc);