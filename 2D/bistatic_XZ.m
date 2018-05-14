%% 基于傅里叶变换反演的双站程序，采用傅里叶变换时一定要满足采样定理，理论上空间谱的网格要均匀划分


close all;
clear all;
clc;

Dx = 0.4;
Dy = 0.4;
objecty_center = 0.34;
objectx_center = 0;
R0 = abs(objecty_center);

%% 发射天线设置，发射机固定，接收机动
c = 3e8;
B = 10e9;
f0 = 180e9;
lambda = c/(f0+B/2);
lambda_min=c/f0;
coeff  = 1;                        % 频域过/降采样系数
fsampling = c/(2*Dy)*coeff;      
fvec = (f0-B/2:fsampling:f0+B/2).';
kvec = 2*pi*fvec/c;
fnum = length(fvec);
Tx = 0.2;
Ty = 0;

%% 接收天线设置
Lx = 0.2;
coefRx = 1;                         % 方位向过/降采样系数
Rxsampling = lambda_min/4*coefRx;              
Ry = 0;
Rx = (-Lx/2:Rxsampling:Lx/2).';
Rxnum = length(Rx);

%% 参考天线设置
refx = Dx + objectx_center;
refy = -Dy/2 + objecty_center;
refamp = 1;

%% 目标场景设置
if Lx < Dx
    theta_span = 2 * asin((Dx/2)/sqrt(R0^2+(Dx/2)^2));
else
    theta_span = 2 * asin((Dx/2+Lx/2)/sqrt(R0^2+(Dx/2+Lx/2)^2));
end

coefbistatic = 2;           % 双站分辨率下降因子
rhox = lambda/4/sin(theta_span/2);
rhoy = c/2/B;
coefx = 6;
coefy = 2;
% object = [-3 -3 1
%           -2 -2 1    
%           -1 -1 1
%           0 0 1
%           1 1 1
%           2 2 1
%           3 3 1
%           -1 1 1
%           -2 2 1
%           -3 3 1
%           1 -1 1
%           2 -2 1
%           3 -3 1
%           1 0 1
%           2 0 1
%           3 0 1
%           -1 0 1
%           -2 0 1
%           -3 0 1
%           0 1 1
%           0 2 1
%           0 3 1
%           0 -1 1
%           0 -2 1
%           0 -3 1] * diag([coefx*rhox coefy*rhoy 1]);
%% 绘制金属平面
X0=0;
Y0=0;
LX=2;
LY=1;
NX=5;
NY=4;
X1=linspace(-LX/2,LX/2,NX);
% Y1=(y0+Ly/2)*ones(1,NX);
Y2=(Y0-LY/2)*ones(1,NX);

Y3=linspace(-LY/2,LY/2,NY);
X2=(X0+LX/2)*ones(1,NY);
X3=(X0-LX/2)*ones(1,NY);
X=[X1,X2,X3];
Y=[Y2,Y3,Y3];
Xa2=X1+3;
Ya2=Y2+2;
OBJECT=[0 0;
%     0 0.5
    0 1
%     0 1.5
    0 2
    1 0 
%     1.5 0
    2 0
%     0 -0.5
    0 -1
%     0 -1.5
    0 -2
    1 2 
%     1.5 2
    2 2
    ];
object = [OBJECT ones(1,length(OBJECT(:,1))).' ;
    ] * diag([coefx*rhox coefy*rhoy 1]);

[tempx,tempx_value] = size(find(abs(object(:,1))>Dx/2));
[tempy,tempy_value] = size(find(abs(object(:,2))>Dy/2));
if (tempx == 0 & tempy==0)<1    
    disp('ERROR：目标超过区间设定');
    return;
end

object(:,1) = object(:,1) + objectx_center;      
object(:,2) = object(:,2) + objecty_center;
objectnum = length(object(:,1));
% objectnum = 1;



%% 回波建立
j=sqrt(-1);
dis_ref = sqrt((Tx-refx)^2+(Ty-refy)^2);
s_ref = refamp * exp(-j*kvec*dis_ref);

S = [];
svec = [];
for index_rx = 1:Rxnum
    s = 0;
    for index_object = 1:objectnum
        x = object(index_object,1);
        y = object(index_object,2);
        amp = object(index_object,3);
        dis_tx = sqrt((x-Tx)^2+(y-Ty)^2);
        dis_rx = sqrt((x-Rx(index_rx))^2+(y-Ry)^2);
        s = s + exp(-j*kvec*(dis_tx+dis_rx));
    end
    s = s.*conj(s_ref);
    S =[S;s.'];
    svec = [svec;s];
end
% figure,plot(Rx,abs(S));
      
%% 波数域建立
% *** 理论上这个波数域是三维空间的一个非均匀分布矩阵，但是为了在后续做傅里叶变换...
% *** 只能提取其中的固定kx和ky，否则无法用目标场景的离散网格去做变维度的逆傅里叶变换...
% *** 这是不是说明如果直接采用IFT公式，做恒定维度的逆傅里叶变换，就可以利用这个非均匀分布矩阵

Kx = [];
Kx_adjust = [];
Ky = [];
Ky_adjust = [];
S_ftx = [];
%     figure
for index_f = 1:fnum
    k = kvec(index_f);
    kxmax = k*sin(theta_span/2);                % kx = 2k;
    kxvec = (linspace(-kxmax,kxmax,Rxnum).');
    s = S(:,index_f);
    s_ftx = exp(-j*kxvec*Rx.')*s;               % 对回波方位向做傅里叶变换
    S_ftx = [S_ftx s_ftx];
    kxvec_adjust = kxvec + k*(objectx_center-Tx)/sqrt((objectx_center-Tx)^2+(objecty_center-Ty)^2);
    kyvec = sqrt(k^2-kxvec.^2);
    kyvec_adjust = kyvec + k*(objecty_center-Ty)/sqrt((objectx_center-Tx)^2+(objecty_center-Ty)^2);
    Kx = [Kx kxvec];
    Kx_adjust = [Kx_adjust kxvec_adjust];
    Ky = [Ky kyvec];
    Ky_adjust = [Ky_adjust kyvec_adjust];

%     plot3(kxvec,kyvec,sqrt(kxvec.^2+kyvec.^2),'*');hold on;
%     plot3(kxvec_adjust,kyvec_adjust,sqrt(kxvec_adjust.^2+kyvec_adjust.^2),'*');hold on;
end

% grid on;xlabel('kx'),ylabel('ky'),zlabel('k');hold off;



%% 回波做相位补偿

Phasecomp = ones(Rxnum,1)* (kvec.')*(sqrt((objectx_center-Tx)^2+(objecty_center-Ty)^2)...
            -(objecty_center-Ty)*objecty_center/sqrt((objectx_center-Tx)^2+(objecty_center-Ty)^2)...
            -(objectx_center-Tx)*objectx_center/sqrt((objectx_center-Tx)^2+(objecty_center-Ty)^2)...    
            -dis_ref) + Ky.* Ry;  

S_comp =  S_ftx.*exp(j*Phasecomp);

%% 对目标反演
% 
%% 2-对相位补偿后的回波做两维维度不变的逆傅里叶变换
% xgrid = (rhox/2*((0:Rxnum-1)-Rxnum/2) + objectx_center).';
% ygrid = (rhoy/2*((0:fnum-1)-fnum/2) + objecty_center).';
% 
% 
% S_ifty = [];
% for index_rx = 1:Rxnum
%     ky_adjust = (Ky_adjust(index_rx,:)).';
%     temp1 = (S_comp(index_rx,:)).';
%     temp2 = exp(j*ygrid*(ky_adjust.'))*temp1;
%     S_ifty = [S_ifty;temp2.'];
% end
% 
% S_iftxy = [];
% for index_f = 1:fnum
%     kx_adjust = Kx_adjust(:,index_f);
%     temp1 = S_ifty(:,index_f);
%     temp2 = exp(j*xgrid*(kx_adjust.'))*temp1;
%     S_iftxy = [S_iftxy temp2];
% end
% figure,imagesc(xgrid,ygrid,abs((S_iftxy).'));
% hold on
% plot(object(:,1),object(:,2),'ro');
% figure
% surf(abs(S_iftxy));
% for index_x = 1:length(xgrid)
%     for index_y = 1:length(ygrid)
%         if(abs(S_iftxy(index_x,index_y))<5.8e+06)
%             S_iftxy(index_x,index_y) = 0;
%         end
%     end
% end
% figure,imagesc(xgrid,ygrid,abs((S_iftxy).'));
% hold on
% plot(object(:,1),object(:,2),'ro');
% figure
% surf(abs(S_iftxy));

%% 3-对相位补偿后的回波自行做两维逆傅里叶变换
xstep = rhox/8;
ystep = rhoy/8;
xgrid = (-Dx/2:xstep:Dx/2).' + objectx_center;
ygrid = (-Dy/2:ystep:Dy/2).' + objecty_center;

kx = Kx_adjust(:,round((fnum-1)/2)+1);
kx_adjust = Kx_adjust(:,round((fnum-1)/2)+1);
ky_adjust = (Ky_adjust(round(Rxnum/2),:)).';
 
S_comp_fixedkx = S_comp;
[XX YY] = meshgrid(xgrid,ygrid);

S_ifty = [];
for index_x = 1:Rxnum
    temp = exp(j*ygrid*(Ky_adjust(index_x,:)))* (S_comp_fixedkx(index_x,:).');
    S_ifty = [S_ifty;temp.'];
end
S_iftxy = [];
for index_f = 1:length(ygrid)
    temp = exp(j*xgrid*(kx_adjust.'))* S_ifty(:,index_f);
    S_iftxy = [S_iftxy temp];
end

figure,
% subplot(2,1,1),surf(XX,YY,(abs((S_iftxy).')));
% title('重构图像（三维）');
% xlabel('x(m)'),ylabel('y(m)'),zlabel('E');
imagesc(xgrid*100,ygrid*100,mat2gray(abs((S_iftxy).'))); 
set(gca, 'YDir', 'normal')
% title('重构图像（二维）');
% xlabel('x(m)'),ylabel('y(m)');
% hold on
% plot(object(:,1)*100,object(:,2)*100,'ro');
% legend('为目标真实位置')
% axis([-5 5 40 60])
hold on 
plot(object(:,1)*100,object(:,2)*100,'ro','LineWidth',1.5,'MarkerSize',6);
xlabel('Azimuth(cm)','Fontname','Times New Roman','FontSize',14)
ylabel('Range(cm)','Fontname','Times New Roman','FontSize',14)
% axis([-5 5 40 60])
set(gca,'Fontname','Times New Roman','FontSize',14)


% for index_i = 1:length(xgrid)
%     for index_j = 1:length(ygrid)
%         if abs(S_iftxy(index_i,index_j))<7e+06
%             S_iftxy(index_i,index_j) = 0;
%         end
%     end
% end
% 
% 
% figure,
% imagesc(xgrid,ygrid,mat2gray(abs((S_iftxy).'))); 
% hold on
% plot(object(:,1),object(:,2),'ro');
% figure
% surf(abs(S_iftxy));

%% 4 对距离维做离散化压缩，再进行方位维离散化
% xstep = rhox;
% ystep = rhoy;
% xgrid = (-Dx/2:xstep:Dx/2).' + objectx_center;
% ygrid = (-Dy/2:ystep:Dy/2).' + objecty_center;
% 
% S_comp_fixedkx = S_comp;
% 
% target_x = [];
% for index_x = 1:Rxnum
%     Fy = exp(-j*Ky_adjust(index_x,:)'*(ygrid'));
%     Fyy = Fy' * Fy;
%     Fyyinv = inv(Fyy);
%     temp = Fyyinv * Fy' * S_comp_fixedkx(index_x,:).';
%     target_x = [target_x;temp.'];
% end
% 
% kx_adjust = Kx_adjust(:,round((fnum-1)/2)+1);
% Fx = exp(-j*kx_adjust*(xgrid.'));
% Fxx = Fx.' * conj(Fx);
% Fxxinv = inv(Fxx);
% target = target_x.' * conj(Fx) *Fxxinv;
% 
% figure,
% % subplot(2,1,1),imagesc(xgrid.',ygrid.',abs(target_x.'));
% imagesc(xgrid.'*100,ygrid.'*100,mat2gray(abs(target)));
% set(gca, 'YDir', 'normal')
% hold on
% plot(object(:,1)*100,object(:,2)*100,'ro');
% axis([-5 5 40 60])
% figure
% surf(abs(target));
% 
% for index_x = 1:length(xgrid)
%     for index_y = 1:length(ygrid)
%         if(abs(target(index_y,index_x))<55)
%             target(index_y,index_x) = 0;
%         end
%     end
% end
% figure,
% imagesc(xgrid,ygrid,mat2gray(abs(target)));
% hold on
% plot(object(:,1),object(:,2),'ko');
% xlabel('x(m)','FontSize',14,'FontName','Times New Roman'),ylabel('y(m)','FontSize',14,'FontName','Times New Roman');
% legend({'目标真实位置'},'FontSize',14);
% colormap(gray);