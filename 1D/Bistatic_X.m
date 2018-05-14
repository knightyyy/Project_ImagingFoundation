%% 一维方向维孔径成像
clc;
clear all;
close all;
%% 已知条件
c=3e8;
fc = 183e9;                       % 雷达发射信号中心频率
B = 0;                       % 带宽
lambda = c/fc;                  % 中心频率波长
RZ = -0.44;                      % 目标中心距离与探测器移动平面的最大距离0.5m
k = 2 * pi * fc/c;               % 波数域频率
% L=0.1;%合成孔径一半的距离10cm
thetaX_annt= 6 * pi/180;        % 天线方向角

%% 合成孔径设置
% D=2*L; %孔径实际大小 
RX_num =256;                                      % 探测器方位维的采样点数，设置为奇数
OverSamplingX = 0.5;                            % OverSamplingAzimuth = 1，Nyquist Sampling；> 1 过采样；  < 1 降采样（此时，注意距离模糊）
deltaX = lambda/2/OverSamplingX;             % 探测器方位维采样间隔
RX_pos = vec((-(RX_num-1)/2:(RX_num-1)/2)*deltaX);  % 探测器方位维的采样坐标，     
Lx = (RX_num-1)*deltaX/2;                         % 方位维一半的孔径 
thetaX_span = 2 * asin(Lx/sqrt(RZ^2+Lx^2));         % 目标方位向距离和探测器方位向展开确定的展开角度，
thetaX = min(thetaX_annt, thetaX_span);                % 确定波数域的展开角度

TX = Lx/2;
TZ = -0.5;


%% 设定目标点所在位置的有效范围
% X方位向分辨率
rho_x = abs(RZ)*lambda/(2*(2*Lx));                   % 目标方位向分辨率
ObjectX = (2*Lx + 2*abs(RZ)*tan(thetaX/2)) *2 * OverSamplingX;    % 目标方位向的最大可测距离
FineX = 4;
ObjectX_num = FineX * floor(ObjectX/rho_x);             % 目标方位向的网格数
% ObjectX_num = Nx; 
ObjectX_pos = vec(linspace(-ObjectX/2,ObjectX/2,ObjectX_num));   %目标区域x范围的网格坐标

% 设定方位向距离
Ptar=[
      0,1
      1,1
      2,1
 ]*diag([3*rho_x 1]);
Object_num = length(Ptar(:,1));
% 

%% 设定参考通道坐标
RefX = ObjectX;
RefZ = 0;

%% 基于指数信号构建直达波和回波
% phi_ref = 2*pi*randn(1);
% phi_echo = 2*pi*randn(1);
% 
% sigma_ref = 1;
% R_ref = sqrt((RefX-TX).^2 +(RefZ-TZ).^2);
% s_ref = sigma_ref*exp(1j*k*R_ref+1j*phi_ref);
% 
% S = zeros(RX_num,1);
% for index_x=1:RX_num;
%     s = 0;    
%     for i=1:1:Object_num;
%         xn = Ptar(i,1);
%         sigma = Ptar(i,2);
%         R_TX = sqrt((xn - TX).^2 +TZ.^2);
%         R_RX = sqrt((xn - RX_pos(index_x)).^2 +RZ.^2);
%         s = s + sigma*exp(1j * k * (R_TX+R_RX)+1j*phi_echo);
%     end
%     s = s.* conj(s_ref);
%    S(index_x,:) = s;
% end
% figure,stem(real(S)),hold on, stem(imag(S),'r'),stem(abs(S),'k');


%% 利用正弦接收回波，
f_if = 1e9;
f_lo = 1.1e9;         % f_lo > f_if, 差频为正,
fs = 10e9;
t0 = randi(10,1)*1/fs;
% t0 = 0;
Tnum = 1000;           % 4的整数倍
t = t0 + (0:Tnum)*1/fs;
phi_ref = 2*pi*randn(1);
phi_echo = 2*pi*randn(1);
% phi_ref = 0;
% phi_echo = 0;


b = fir1(48,[0.001 0.2]);
% figure,freqz(b,1,512)

sigma_ref = 1;
R_ref = sqrt((RefX-TX).^2 +(RefZ-TZ).^2);
s_ref_real = sigma_ref * cos(2*pi*f_if*t - k * R_ref - phi_ref);
s_ref_real_noise = awgn(s_ref_real,1);
% 正交变换
s_ref_in = s_ref_real_noise.* cos(2*pi*f_lo*(t-t0));
s_ref_quad = s_ref_real_noise.* sin(2*pi*f_lo*(t-t0));


shape = 'same'; % full,same,valid
s_ref_in = conv(s_ref_in,b,shape);
s_ref_quad = conv(s_ref_quad,b,shape);

s_ref = s_ref_in + 1j * s_ref_quad;
% figure, subplot(2,1,1),plot(s_ref_real),hold on,plot(s_ref_in,'r');
%                        hold on,plot(s_ref_quad,'k');
%                        legend('DirectCos','In','Quad')
%         subplot(2,1,2),plot(abs(fft(s_ref_real))),hold on,plot(abs(fft(s_ref_in)),'r');
%                        hold on,plot(abs(fft(s_ref_quad)),'y');
%                        hold on,plot(abs(fft(s_ref)),'k');

S = zeros(RX_num,1);
Scope = Tnum/2;         % 选择取平均的范围


for index_x=1:RX_num;
    s = 0;    
    for i=1:1:Object_num;
        xn = Ptar(i,1);
        sigma = Ptar(i,2);
        R_TX = sqrt((xn - TX).^2 +TZ.^2);
        R_RX = sqrt((xn - RX_pos(index_x)).^2 + RZ.^2);
        s = s + sigma*cos( 2*pi*f_if*t - k * (R_TX+R_RX) - phi_echo);
%         s = s + sigma*exp(1j*2*pi*(f_if-f_lo)*t+1j * k * (R_TX+R_RX)+1j*phi_echo);
    end
    s_noise = awgn(s,1);
    s_in = s_noise.* cos(2*pi*f_lo*(t-t0));
    s_quad = s_noise.* sin(2*pi*f_lo*(t-t0));
    s_in = conv(s_in,b,shape);
    s_quad = conv(s_quad,b,shape);
    s_exp = s_in + 1j*s_quad;
    s_downsampling = s_exp.* conj(s_ref);                  
    % ******* 特别需要注意这个地方的符号，要保证差频后 s_exp 和 s_ref 的相位符号为﹢************** %
    s_mean = mean(s_downsampling(Tnum/2-Scope/2:Tnum/2+Scope/2));           
    S(index_x,:) = s_mean;
end
% figure,plot(real(s_downsampling),'o'),hold on,plot(imag(s_downsampling),'ro');
% % figure,
% subplot(2,1,1),plot(s),               
%                hold on,plot(s_in,'r');
%                hold on,plot(s_quad,'k');
%                legend('EchoCos','In','Quad');
%                hold on, plot(s_noise,'g');
% % subplot(2,1,2),plot(abs(fft(s))),hold on,plot(abs(fft(s_in)),'r');
% %                hold on,plot(abs(fft(s_quad)),'y');
% %                hold on,plot(abs(fft(s_exp)),'k');
% % figure,stem(real(S)),hold on, stem(imag(S),'r'),stem(abs(S),'k');

SNR = 5;
S = awgn(S,SNR);

% 波数域频率
kx_range = k * Lx /(sqrt((Lx).^2+RZ.^2));
kx = vec(linspace(-kx_range,kx_range,RX_num));
kx_adjust = kx + k.*TX/sqrt(TX^2+TZ^2);

kz = sqrt(k.^2 - kx.^2);
kz_adjust = kz*abs(RZ) + k*sqrt(TX^2+TZ^2) - k*abs(R_ref);


%% 对回波进行一维傅里叶变换
% 对方位向x进行傅里叶变换
S_ftx = exp(-1j * kx * RX_pos.') * S;

% 进行相位补偿
S_comp = S_ftx.* exp(-1j * kz_adjust);
figure
plot(abs(S_comp));
% 对X，Y方位向进行傅里叶逆变换
%对方位向x进行傅里叶逆变换
S_iftx = exp(1j * ObjectX_pos * kx_adjust.') * S_comp./RX_num;
S_iftx=S_iftx./(max(S_iftx));

%% 直接用IFT
figure
% subplot(1,2,1),surf(XX,YY,abs(S_iftxy).'),title('2D Sum');
% hold on, plot3(Ptar(:,1),Ptar(:,2),Ptar(:,3),'wo','MarkerSize',6)
plot(ObjectX_pos*500,abs(S_iftx)),xlabel('azimuth/cm'),ylabel('amplitude');
hold on, 
% plot(Ptar(:,1)*500,1,'ko','LineWidth',1.5,'MarkerSize',5);
stem(Ptar(:,1)*500,Ptar(:,2),'fill','-.')
% stem(Ptar(2,1)*500,Ptar(2,2),'fill','-.')
% set(gca,'XTickMode','manual','XTick',[Ptar(1,1),Ptar(2,1),Ptar(3,1)]);
axis([-15 15 0 1]);
% grid on
% title()
legend('Measured','True target')


