%% 柱面成像算法0103版(完成版)
clc;
clear all;
close all;
%% 目标图片读取
% 
% A=imread('C:\Users\dell\Desktop\中文小论文图片\2.jpg ');
% t=graythresh(A); 
% Y=im2bw(A,t);
% q=8;
% Y1=Y(1:q:end,1:q:end);
% figure
% imshow(Y1);
% [Ptar_x1,Ptar_y1]=find(Y1==0);
% B=imread('C:\Users\dell\Desktop\中文小论文图片\4.jpg ');
% t=graythresh(B); 
% Y=im2bw(B,t);
% Y2=Y(1:q:end,1:q:end);
% figure
% imshow(Y2);
% [Ptar_x2,Ptar_y2]=find(Y2==0);
%% 参数设置
c=3e8;%光速
fc=35e9;%发射信号频点
lambda=c/fc;%中心频率波长
k=2*pi*fc/c;%波数域频率
R_det=1;%观测柱面到圆心的距离
length_z=2;%观测平面的高度
R_tr0=0.2;%目标柱面到达圆心的距离
%% 角度维采样间隔计算
N_theta=720;%角度维采样个数
theta_range=pi;
delta_theta=theta_range/(N_theta);%角度维采样间隔
theta_deg=delta_theta*180/pi;
theta_det=(-N_theta/2:(N_theta/2-1))*delta_theta;%水平方向上探测器的扫描角度
x_det=R_det*cos(theta_det);%采样点的横坐标位置
y_det=R_det*sin(theta_det);%采样点的纵坐标位置
%% 高度维采样间隔计算
% delta_z=lambda*sqrt((R_det^2)+(length_z/2)^2)/(length_z);%竖直平面采样位置
delta_z=lambda/3;
N_z=floor(length_z/delta_z);
z_det=(0:(N_z-1))*delta_z;%竖直方向上探测器采样位置
%% 目标点参数设置
Ptar=[
50 100 R_tr0
0 100 R_tr0
% -50 100 R_tr0
-50 50 R_tr0
-50 0 R_tr0
50 0 R_tr0
% -50 150 R_tr0
-50 200 R_tr0
0 200 R_tr0
50 200 R_tr0
0 0 R_tr0
 ]*diag([2*delta_theta 3*delta_z 1]);
% Ptar=[
% % [Ptar_x1.'-80 Ptar_x2.'].'  [Ptar_y1.' Ptar_y2.'].' [(R_tr0)*ones(1,(length(Ptar_x1))) (R_tr0)*ones(1,(length(Ptar_x2)))].'
% ]*diag([5*delta_theta 5*delta_z 1]);
Ntar=size(Ptar,1);
%% 回波信号构建，探测器同一位置自发自收
for index_det=1:N_theta;
    s=0;
    for index_tr=1:Ntar;
        theta_tr=Ptar(index_tr,1);
        z_tr=Ptar(index_tr,2);
        R_tr=Ptar(index_tr,3);
        sigma=1;%目标点的散射系数
        R=sqrt((R_tr*cos(theta_tr)-x_det(index_det)).^2+(R_tr*sin(theta_tr)-y_det(index_det)).^2+(z_tr-z_det).^2);%收发机与目标之间的直线距离
        phase_diff=2*k*R;%相位
        s=s+sigma*exp(-1j*phase_diff);%回波信号
    end
    ECHO(index_det,:)=s;
end
% 
% figure
% % plot(x_det,y_det,'*');
%  imagesc(abs(ECHO));
%% 设定两个中间变量
Ptar1=[0 0 R_tr0];
index_tr=1;
% for index_tr=1:Ntar;
m=(R_det-Ptar1(index_tr,3)).^2+(Ptar1(index_tr,2)-z_det).^2;%1*N
n=R_det*Ptar1(index_tr,3);
q=(R_det-Ptar1(index_tr,3)).^2;
phi=Ptar1(index_tr,1)-theta_det;

%% 确定波数域频率(非均匀)
for index_z=1:N_z;
    for index_phi=1:N_theta;
        ktheta_origin(index_phi,index_z)=2*k*n*sin(phi(index_phi))/(m(index_z)+n*(phi(index_phi))^2);
        kh_origin(index_phi,index_z)=2*k*(Ptar1(index_tr,2)-z_det(index_z))/(m(index_z)+n*(phi(index_phi))^2);
    end
end
%% 确定波数域频率（均匀）由于每一个圆周上的结果都是单调递增(均匀网格1)
for index_z=1:N_z;
    ktheta(:,index_z)=linspace(min(min(ktheta_origin)),max(max(ktheta_origin)),N_theta);
end
for index_theta=1:N_theta;
    kh(index_theta,:)=linspace(min(min(kh_origin)),max(max(kh_origin)),N_z);
end
%% 确定波数域频率（均匀）由于每一个圆周上的结果都是单调递增(均匀网格2)
% for index_z=1:N_z;
%     ktheta(:,index_z)=linspace(min(ktheta(:,index_z)),max(ktheta(:,index_z)),N_theta);
% end
% for index_theta=1:N_theta;
%     kh(index_theta,:)=linspace(min(kh(index_theta,:)),max(kh(index_theta,:)),N_z);
% end
%% 二维傅里叶变换
for index_z=1:N_z
    s_FTx=exp(-1j*(ktheta(:,index_z)*theta_det))*ECHO(:,index_z);
    S_FTy(:,index_z)=s_FTx;
end
for index_x=1:N_theta;
    s_FTz=exp(-1j*(kh(index_x,:)).'*z_det)*(S_FTy(index_x,:).');
    S_FTz(index_x,:)=s_FTz.';
end
% figure
% imagesc(abs(S_FTz));
%% 进行相位补偿
for index_z=1:N_z;
    for index_phi=1:N_theta;
        k1(index_phi,index_z)=sqrt(4*k^2-(ktheta(index_phi,index_z)^2)/n);
        phase(index_phi,index_z)=sqrt(k1(index_phi,index_z)^2-(kh(index_phi,index_z))^2)*sqrt(q);
    end
end
 S_comp=S_FTz.*exp(1j*phase);
%  figure
%  sigma_theta=pi/max(max(ktheta));
% theta_grid=linspace(-theta_range/2,theta_range/2,N_theta);
% sigma_z=lambda/2;
% z_grid=linspace(-1,length_z,N_z);
% imagesc(180*theta_grid/pi,z_grid,abs(S_comp).');
% hold on
% plot(Ptar(:,1)*180/pi,Ptar(:,2),'ro');
%% 采用二阶巴特沃斯低通滤波器
[M_comp,N_comp]=size(S_comp);
nn=2;%二阶
d0=50;
m_comp=fix(M_comp/2);
n_comp=fix(N_comp/2);
for index_M=1:M_comp;
    for index_N=1:N_comp;
        d_comp=sqrt((index_M-m_comp)^2+(index_N-n_comp)^2);
        if d_comp==0;
            h=0;
        else
            h_comp=1/(1+(d_comp/d0)^(2*nn));
%             h_comp=exp(-(d_comp*d_comp)/(2*d0*d0));
        end
        S_comp(index_M,index_N)=h_comp*S_comp(index_M,index_N);
    end
end
%  figure
%  sigma_theta=pi/max(max(ktheta));
% theta_grid=linspace(-theta_range/2,theta_range/2,N_theta);
% sigma_z=lambda/2;
% z_grid=linspace(-1,length_z,N_z);
% imagesc(180*theta_grid/pi,z_grid,abs(S_comp).');
% hold on
% plot(Ptar(:,1)*180/pi,Ptar(:,2),'ro');
% %% 进行stolt插值处理
% DKtheta=abs(2*k-sqrt(4*k^2-(ktheta.^2)/n));
% P=8/2;%8点sinc插值
% NDKtheta=fix(DKtheta);
% Nd=NDKtheta(N_theta,1);
% B=[S_comp S_comp(:,1:Nd+P)];
% S_comp1=zeros(N_theta,N_z);
% for index_phi=1:N_theta;
%     for index_z=P:N_z;
%         NN=NDKtheta(index_phi,index_z)+(-P+1:P);
%         S_comp1(index_phi,index_z)=sinc(DKtheta(index_phi,index_z)-NN)*B(index_phi,index_z+NN).';
%     end
% end
%% 绘制成像结果网格(二维傅里叶逆变换)
%% planA
% sigma_theta=pi/max(max(ktheta));
% theta_grid=linspace(-pi,pi,N_theta);
% sigma_z=lambda/2;
% N_z_grid=length_z/sigma_z;
% z_grid=linspace(0,length_z,N_z_grid);
% %% 进行二维傅里叶逆变换
% for index_theta=1:N_theta;
%     s_IFTz(index_theta,:)=(exp(1j*(z_grid).'*kh(index_theta,:))*(S_comp(index_theta,:).'))./N_theta;
% end
% for index_theta=1:N_theta;
%     ktheta_new(index_theta,:)=linspace(min(ktheta(index_theta,:)),max(ktheta(index_theta,:)),N_z_grid);
% end
% for index_z=1:N_z_grid;
%     S_IFTx(index_tr,:,index_z)=(exp(1j*(2*theta_grid).'*ktheta_new(:,index_z).')*(s_IFTz(:,index_z)))./N_z_grid;
% end
%% planB
sigma_theta=pi/max(max(ktheta));
theta_grid=linspace(-theta_range/2,theta_range/2,N_theta);
sigma_z=lambda/2;
z_grid=linspace(-1,length_z,N_z);
%% 进行二维傅里叶逆变换
for index_z=1:N_z;
    S_IFTx(:,index_z)=(exp(1j*(theta_grid).'*ktheta(:,index_z).')*(S_comp(:,index_z)))./N_z;
end
figure
imagesc(180*theta_grid/pi,z_grid,abs(S_IFTx).');
hold on
plot(Ptar(:,1)*180/pi,Ptar(:,2),'ro');
for index_theta=1:N_theta;
    s_IFTz(index_theta,:)=(exp(1j*(z_grid).'*kh(index_theta,:))*(S_IFTx(index_theta,:).'))./N_theta;
end
figure
imagesc(180*theta_grid/pi,z_grid,abs(s_IFTz).');
set(gca, 'YDir', 'normal');
hold on
plot(Ptar(:,1)*180/pi,Ptar(:,2),'ro','LineWidth',2,'MarkerSize',10);
set(gca, 'YDir', 'normal');
xlabel('Degree(°)','Fontname','Times New Roman','FontSize',14);
ylabel('Height(m)','Fontname','Times New Roman','FontSize',14);
axis on
set(gca,'Fontname','Times New Roman','FontSize',14)
%% planC 采用压缩感知的方式进行处理
% A=[];
% for index_Rz=1:N_z;
%     for index_Rx=1:N_theta;
%         ktheta_CS=ktheta(index_Rx,index_Rz);
%         kh_CS=kh(index_Rx,index_Rz);
%         for index_z=1:N_z;
%             Anm=[];
%             z=z_grid(index_z);
%             for index_theta=1:N_theta;
%                 theta=theta_grid(index_theta);
%                 nm=(index_theta-1)*length(z_grid)+index_z;
%                 Anm(nm)=exp(-1j*ktheta_CS*theta-1j*kh_CS*z);
%             end
%         end
%         A=[A;Anm];
%     end
% end
% end

% for index_tr=1:Ntar;
% figure
%     S_change(:,:)=S_IFTx(index_tr,:,:);
% %     figure(index_tr);
%     imagesc(180*theta_grid/pi,z_grid,abs(S_change).');
%     hold on
%     plot(Ptar(:,1)*180/pi,Ptar(:,2),'ro');
%     phase = angle(S_change);
% %     figure
% %      imagesc(180*theta_grid/pi,z_grid,angle(S_change).');
% %      hold on,  plot(Ptar(:,1)*180/pi,Ptar(:,2),'ro');
% %      figure,surf(phase)
% % end
