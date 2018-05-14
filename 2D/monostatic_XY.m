%% ��ά������ά����ά
clc;
clear all;
close all;
%% ��֪����
c=3e8;
fc=180e9;                       % �״﷢���ź�����Ƶ��
B = 0;                       % ����
lambda = c/fc;                  % ����Ƶ�ʲ���
TRZ = - 0.5;                      % Ŀ�����ľ�����̽�����ƶ�ƽ���������0.5m
k = 2 * pi * fc/c;               % ������Ƶ��
% L=0.1;%�ϳɿ׾�һ��ľ���10cm
thetaX_annt= 20 * pi/180;        % ���߷����
thetaY_annt= 20 * pi/180;        % ���߷����
%% �ϳɿ׾�����
% D=2*L; %�׾�ʵ�ʴ�С 
TRX_num = 257;                                      % ̽������λά�Ĳ�������������Ϊ����
OverSamplingX = 1;                            % OverSamplingAzimuth = 1��Nyquist Sampling��< 1 ��������   > 1 ����������ʱ��ע�����ģ����
deltaX = lambda/2/OverSamplingX;             % ̽������λά�������
TRX_pos = vec((-(TRX_num-1)/2:(TRX_num-1)/2)*deltaX);  % ̽������λά�Ĳ������꣬     
Lx = (TRX_num-1)*deltaX/2;                         % ��λάһ��Ŀ׾� 
thetaX_span = 2 * asin(Lx/sqrt(TRZ^2+Lx^2));         % Ŀ�귽λ������̽������λ��չ��ȷ����չ���Ƕȣ�
thetaX = min(thetaX_annt, thetaX_span);                % ȷ���������չ���Ƕ�

TRY_num = 257;                                      % ̽������λά�Ĳ�������������Ϊ����
OverSamplingY = 1;
deltaY = lambda/2/OverSamplingY;             % ̽������λά�������
TRY_pos = vec((-(TRY_num-1)/2:(TRY_num-1)/2)*deltaY);  % ̽������λά�Ĳ������꣬     
Ly = (TRY_num-1)*deltaY/2;                         % ��λάһ��Ŀ׾� 
thetaY_span = 2 * asin(Ly/sqrt(TRZ^2+Ly^2));         % Ŀ�귽λ������̽������λ��չ��ȷ����չ���Ƕȣ�
thetaY = min(thetaY_annt, thetaY_span);                % ȷ���������չ���Ƕ�

%% �趨Ŀ�������λ�õ���Ч��Χ
% X��λ��ֱ���
rho_x = abs(TRZ)*lambda/(2*(2*Lx));                   % Ŀ�귽λ��ֱ���
ObjectX = (2*Lx + 2*abs(TRZ)*tan(thetaX/2)) / OverSamplingX;    % Ŀ�귽λ������ɲ����
FineX = 4;
ObjectX_num = FineX * floor(ObjectX/rho_x);             % Ŀ�귽λ���������
% ObjectX_num = Nx; 
ObjectX_pos = vec(linspace(-ObjectX/2,ObjectX/2,ObjectX_num));   %Ŀ������x��Χ����������

rho_y = abs(TRZ)*lambda/(2*(2*Ly));                   % Ŀ�귽λ��ֱ���
ObjectY = (2*Ly + 2*abs(TRZ)*tan(thetaY/2)) / OverSamplingY;    % Ŀ�귽λ������ɲ����
FineY = 4;
ObjectY_num = FineY * floor(ObjectY/rho_y);             % Ŀ�귽λ���������
% ObjectX_num = Nx; 
ObjectY_pos = vec(linspace(-ObjectY/2,ObjectY/2,ObjectY_num));   %Ŀ������x��Χ����������
%% �趨��λ�����
data=xlsread('��Ŀ��.xlsx');
i=1;
data_x=data(i,1);
data_y=data(i,2);
data_sigma=data(i,3);
Ptar1=data;
% Ptar1=[0,-4,1
%       0,-3,1
%       0,-2,1
%       0,-1,1
%       0,0,1
%       0,1,1
%       0,2,1
%       0,3,1
%       0,4,1
%       2,4,1
%       4,4,1
%       6,4,1
%       8,4,1
%       2,0,1
%       4,0,1
%       6,0,1
%       -4-4,4,1
%       -3-4,3,1
%       -2-4,2,1
%       -1-4 1 1
%       1-4,-1,1
%       2-4,-2,1
%       -5-4,3,1
%       -6-4,2,1
%       -7-4,1,1
%       -8-4,0,1
%       -9-4,-1,1
%       -10-4,-2,1
%       -2-4,0,1
%       -4-4,0,1
%       -6-4,0,1
%       -8-4,0,1
%       -2,0,1
%       -4,0,1
%       ];
Ptar=Ptar1* diag([4*rho_x 4*rho_y 1]);
Ptar(:,1) = Ptar(:,1) + 5 * rho_x;
Object_num = length(Ptar(:,1));
% PtarZ_min = -0.0187; PtarZ_max = 0.0187;
PtarZ_min = -0.0766; PtarZ_max = 0.0766;
PtarZ = PtarZ_min + (PtarZ_max - PtarZ_min)*rand(Object_num,1);
% PtarZ = 0.0187*ones(Object_num,1);
% PtarZ = zeros(Object_num,1);
% PtarZ(1) = 0.01;
%% 
figure
scatter3(PtarZ,Ptar(:,1),Ptar(:,2),'filled')
ylabel('��λ��x/m');
zlabel('��λ��y/m');
xlabel('������z/m');
title('Ŀ����ƽ��Ķ�ά��');
% for i=1:Object_num;
%     text(Ptar(i,1),Ptar(i,2),PtarZ(i),['(',num2str(Ptar(i,1)),',',num2str(Ptar(i,2)),',',num2str(PtarZ(i)),')']);
% end
% Object_num = 5;
% ObjectX_seq = randperm(ObjectX_num);           % the x-coordination of target
% Ptar_x = vec((ObjectX_seq(1:Object_num) * rho_x/FineX - ObjectX/2)); % the x-coordination of target, column      
% ObjectZ_seq =  randperm(ObjectZ_num);
% Ptar_z = vec((ObjectZ_seq(1:Object_num) * rho_z/FineZ - ObjectZ/2));  % the z-coordination of target, column        
% Ptar_amp = ones(Object_num,1);
% Ptar = [Ptar_x Ptar_z Ptar_amp];



%% �����ز�

S = zeros(TRX_num,TRY_num);
for index_x=1:TRX_num;
    s = 0;
        for i=1:1:Object_num;
            xn = Ptar(i,1);
            yn = Ptar(i,2);
            zn = PtarZ(i);  
            sigma = Ptar(i,3);
            R = sqrt((xn - TRX_pos(index_x)).^2 + (yn - TRY_pos).^2 + (zn-TRZ).^2);
            s = s + sigma*exp(-1j * 2 * k * R);
        end
   S(index_x,:) = s.';
end
figure
plot(abs(s));
SNR = 10;
S = awgn(S,SNR);

% ������Ƶ��
kx_range = k * 2 * Lx /(sqrt((Lx).^2+TRZ.^2));
ky_range = k * 2 * Ly /(sqrt((Ly).^2+TRZ.^2));
kx = vec(linspace(-kx_range,kx_range,TRX_num));
ky = vec(linspace(-ky_range,ky_range,TRY_num));
Kx = kx * ones(1,TRY_num);
Ky = ones(TRX_num,1) * (ky.'); 
Kz = sqrt((2*k).^2 - Kx.^2 - Ky.^2);

delta_kx = 2*kx_range/TRX_num;
delta_ky = 2*ky_range/TRY_num;
Kz_approximation = 2*k - delta_kx^2/4/k*(vec([-(TRX_num-1)/2:(TRX_num-1)/2].^2))*ones(1,TRY_num) - delta_ky^2/4/k*ones(TRX_num,1)*([-(TRY_num-1)/2:(TRY_num-1)/2].^2);

Kz = Kz_approximation;

Kz_max = max(max(Kz));
Kz_min = min(min(Kz));
Z_threshold = pi/(Kz_max-Kz_min);



%% �Իز����ж�ά����Ҷ�任
% �Է�λ��x���и���Ҷ�任
for index_y = 1:TRY_num;
    temp = exp(-1j * kx * TRX_pos.') * S(:,index_y);
    S_ftx(:,index_y) = temp;
end

for index_x = 1:TRX_num;
    temp = exp(-1j * ky * TRY_pos.') * (S_ftx(index_x,:).');
    S_ftx(index_x,:) = temp.';
end

% ������λ����
S_comp = S_ftx.* exp(1j * Kz * abs(TRZ));
% ��X��Y��λ����и���Ҷ��任
%�Է�λ��x���и���Ҷ��任
S_iftx = exp(1j * ObjectX_pos * kx.') * S_comp./TRX_num;
S_iftxy = (exp(1j * ObjectY_pos * ky.') * (S_iftx.')./TRY_num).';

%% ֱ����IFT
[XX,YY] = meshgrid(ObjectX_pos,ObjectY_pos);
% figure
% % subplot(1,2,1),surf(XX,YY,abs(S_iftxy).'),title('2D Sum');
% % hold on, plot3(Ptar(:,1),Ptar(:,2),Ptar(:,3),'wo','MarkerSize',6)
% subplot(1,2,1),imagesc(ObjectX_pos*1000,ObjectY_pos * 1000,mat2gray(abs((S_iftxy).'))),axis xy;
% subplot(1,2,2),imagesc(ObjectX_pos*1000,ObjectY_pos * 1000,mat2gray(abs((S_iftxy).'))),axis xy;
% xlabel('x(mm)'),ylabel('y(mm)');
% hold on, plot(Ptar(:,1)*1000,Ptar(:,2)*1000,'ko','LineWidth',1.5,'MarkerSize',5);

%% ���Z�Ŀ��ܷ�Χ�����ж��������������ؼ����ӣ�������Ӧ�ļ�Ȩϵ��
ZL = -0.1;   % ZL=-ZR=0.1,Znum=16,17
ZR = 0.1;

ZL = PtarZ_min;
ZR = PtarZ_max;
Znum = 33;
TRZ_interval = abs(ZR-ZL)/Znum
TRZ_profile = vec(linspace(ZL,ZR,Znum)+TRZ);

TRZ_interval = 0.0125;
TRZ_profile = vec([ZL:TRZ_interval:ZR])+TRZ;
Znum = length(TRZ_profile);


S_iftxy_profile = [];
for index_z = 1:Znum
    TRZ_temp = TRZ_profile(index_z);
    S_comp = S_ftx.* exp(1j * Kz * abs(TRZ_temp));
    temp1 = exp(1j * ObjectX_pos * kx.') * S_comp./TRX_num;
    temp2 = (exp(1j * ObjectY_pos * ky.') * (temp1.')./TRY_num).';
    temp3 = vec(temp2);
    S_iftxy_profile = [S_iftxy_profile temp3];
%     figure(3)
%     subplot(1,2,1),imagesc(ObjectX_pos*1000,ObjectY_pos * 1000,mat2gray(abs((S_iftxy).'))),axis xy,title('IFT');
%     % hold on, plot3(Ptar(:,1),Ptar(:,2),Ptar(:,3),'wo','MarkerSize',6)
%     subplot(1,2,2),imagesc(ObjectX_pos*1000,ObjectY_pos * 1000,mat2gray(abs((temp2).'))),axis xy,title('Weighted IFT');
%     xlabel('x(mm)'),ylabel('y(mm)');
%     pause(1)
%     close 
end
coeff = ones(Znum,1)./Znum;
S_iftxy_weighted = S_iftxy_profile * coeff;
S_iftxy_weighted = reshape(S_iftxy_weighted,ObjectX_num,ObjectY_num);

figure
subplot(1,2,1),imagesc(ObjectX_pos*1000,ObjectY_pos * 1000,mat2gray(abs((S_iftxy).'))),axis xy,title('δ��Ȩƽ��������ع�ͼ��');
xlabel('x(mm)'),ylabel('y(mm)');
% hold on, plot3(Ptar(:,1),Ptar(:,2),Ptar(:,3),'wo','MarkerSize',6)
subplot(1,2,2),imagesc(ObjectX_pos*1000,ObjectY_pos * 1000,mat2gray(abs((S_iftxy_weighted).'))),axis xy,title('��Ȩƽ��������ع�ͼ��');
xlabel('x(mm)'),ylabel('y(mm)');

%% ����entropy���趨��Ȩϵ��
% S_iftxy_profile_normal = S_iftxy_profile * (diag(1./max(abs(S_iftxy_profile))));
% S_iftxy_prob = abs(S_iftxy_profile_normal)*(diag(1./sum(abs(S_iftxy_profile_normal))));
% S_iftxy_entropy = vec(sum(-S_iftxy_prob.*log(S_iftxy_prob)));
% [Entropy_sort Entropy_index] = sort(S_iftxy_entropy,'descend');
% coeff_entropy = S_iftxy_entropy(Entropy_index)/sum(S_iftxy_entropy);
% S_iftxy_profile_sort = S_iftxy_profile(:,Entropy_index);
% S_iftxy_weighted_entropy = S_iftxy_profile_sort * flipud(coeff_entropy);
% 
% S_iftxy_weighted_entropy = reshape(S_iftxy_weighted_entropy,ObjectX_num,ObjectY_num);
% figure
% subplot(1,3,1),imagesc(ObjectX_pos*1000,ObjectY_pos * 1000,mat2gray(abs((S_iftxy).'))),axis xy,title('IFT');
% xlabel('x(mm)'),ylabel('y(mm)');
% hold on, plot(Ptar(:,1)*1000,Ptar(:,2)*1000,'ko','LineWidth',1.5,'MarkerSize',5);
% 
% subplot(1,3,2),imagesc(ObjectX_pos*1000,ObjectY_pos * 1000,mat2gray(abs((S_iftxy_weighted).'))),axis xy,title('weighted IFT');
% xlabel('x(mm)'),ylabel('y(mm)');
% hold on, plot(Ptar(:,1)*1000,Ptar(:,2)*1000,'ko','LineWidth',1.5,'MarkerSize',5);
% 
% % hold on, plot3(Ptar(:,1),Ptar(:,2),Ptar(:,3),'wo','MarkerSize',6)
% subplot(1,3,3),imagesc(ObjectX_pos*1000,ObjectY_pos * 1000,mat2gray(abs((S_iftxy_weighted_entropy).'))),axis xy,title('Weighted Entropy IFT');
% xlabel('x(mm)'),ylabel('y(mm)');
% hold on, plot(Ptar(:,1)*1000,Ptar(:,2)*1000,'ko','LineWidth',1.5,'MarkerSize',5);

 
