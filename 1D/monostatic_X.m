%% һά����ά�����㷨
clc;
clear all;
close all;
%% ��֪����
c=3e8;%����
fc=183e9;%�״﷢���ź�����Ƶ��
Yc=0.5;%Ŀ����̽�����ƶ�ƽ��ľ���0.5m
v=0.001;%̽�����ƶ��ٶ�1mm/s
L=0.1;%�ϳɿ׾�һ��ľ���10cm
theta=10*pi/180;%���߷����
%% �ϳɿ׾��Ĳ����趨
lambda=c/fc;%����
w=2*pi*fc;
k=2*w/c;%������Ƶ��
D=2*L;%�׾�ʵ�ʴ�С
%% �趨Ŀ�������λ�õ���Ч��Χ
X_range=2*L+2*Yc*tan(theta/2);
%% ��λ��ֱ���
delta_x=Yc*lambda/(2*D);
% Na=2*L/delta_x;%������ĸ���
Na=256;
x=linspace(-L,L,Na);%̽�����˶���Χ
u=x/v;
% du=2*L/v/Na;
% fu=linspace(-1/2/du,1/2/du,Na);%fu������
%% Ŀ������Χ
xa=linspace(-X_range/2,X_range/2,Na);
%% �趨��λ����뼰����ϵ��
% Ntar=2;
% Ptar=[0,1
%     4*delta_x,1
%     10*delta_x,1
%     20*delta_x,1
%     40*delta_x,1
%     ];
Ntar=2;
G=X_range/2/delta_x;
Ptar=[G*rand(1,Ntar)*delta_x].';%�������Ŀ����λ��
% Ntar=20;
% Ptar=linspace(-10*delta_x,10*delta_x,Ntar);
%% �����ز��ź�
s=0;
for i=1:1:Ntar; 
%     xn=Ptar(i,1);
%     sigma=Ptar(i,2);
    xn=Ptar(i);
    sigma=1;
    R=sqrt((Yc^2)+((xn-x).^2));
    s=s+sigma*exp(-1j*k*R);%�ز��ź� ���Ŀ���ƶ���Χ����L,С��X_range/2ʱ����-L,L�ķ�Χ���޷���ʾ����
end
% figure
% plot(abs(fft(s)));
% SNR_db=10; 
% SNR=10^(SNR_db/10); 
% %%%%%%% %�����ź����� 
% xt=s;
% sigpower=0; 
% for i=1:length(xt) 
%     sigpower=sigpower+xt(i).^2;
% end
% %�������� 
% noisepower=sigpower/SNR;%������������ 
% noise=randn(1,length(xt));%������ֵΪ�㣬����Ϊһ�������˹���� 
% noise=noise-repmat(mean(mean(noise)),1,length(xt));%�þ�ֵ���ӽ�����%����У�� 
% AWGNpower=0; 
% for i=1:length(xt) 
%     AWGNpower=AWGNpower+noise(1,length(xt)).^2;%��������������
% end
% K=noisepower/AWGNpower;%��������������ʵ���������������ı�ֵ 
% K1=sqrt(K);%У������������Ȩֵ 
% snr=10*log10(sigpower/sum((K1*noise).^2));%ʵ����������
% s=xt+(K1*noise);
%% �Իز�����ʱ������������Ҷ�任
% Kx=(k*((x-xa))./sqrt((Yc^2)+(((x-xa)).^2)));
kx_range=2*pi*X_range/(lambda*sqrt((X_range/2).^2+Yc.^2));
Kx=linspace(-kx_range,kx_range,Na);
Sn=exp(-1j*Kx.'*x)*s.';
Sn=Sn.';
% Sn=fftshift(fft(fftshift(s)));
%% ������λ����
%     Kx=(k*(xn-x)./sqrt((Yc^2)+((xn-x).^2)));
val=exp(1j*((sqrt(((k.^2))-(Kx.^2))))*Yc);                                                                                                                                                                                                                                                                                                                                                                                            
S_comp=Sn.*val;
%% �Է���λ���з�����Ҷ�任                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
% s_ifft=fftshift(ifft(fftshift(S_comp)));
s_ifft=exp(1j*xa.'*Kx)*S_comp.';
s_ifft=s_ifft.';
%% �Է�ֵ���й�һ������
s_ifft=s_ifft./max(s_ifft);
% %% ��ͼ
% figure                                      
% plot(xa,abs(s_ifft));
% xlabel('λ��/m');
% ylabel('��ֵ');
% axis on;
% s_fft=fftshift(fft(fftshift(s_ifft)));
figure
plot(xa,abs(s_ifft));
hold on
plot(Ptar(:),ones(1,length(Ptar)),'r*');
% title('Ƶ��ͼ');
% vm=sqrt(4-((D.^2)/((L.^2)+Yc.^2)));
% delta_z=c/(2*(2-vm)*fc);