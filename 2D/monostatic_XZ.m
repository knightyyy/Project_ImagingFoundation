%%&%%%%% NUFFT test in thesis P20 %%%%%%%%%%%%%%
%%%% programming on equation
close all;
clc;
clear all;

% capital form: matrix
% capital +lower case: constant
% lower case: variable,vector

%% echo reconstruction

fs=75e9;                % starting frequency
B=35e9;                  % Bandwidth
Nf=80;                  % sampling number in frequency domain
f_step=B/Nf;            % sampling interval in frequency domain
f=fs+(0:Nf-1)*f_step;   % the sequence of frequency
c=3e8;                  % the velocity of light
lambda_min=c/fs;

Na=256;                 % number of antenna array
x_step=lambda_min/2;    % sampling interval in x direction
R0=2;                   % the distance between the array and OBJECT
y_tr=-R0.*ones(Na,1);   % the y-coordination of TR, including sign
x_tr=((0:Na-1).'-Na/2)*x_step;  % the x-coordination of TR, including sign
Lx=(Na-1)*x_step;       % efficient aperture length

theta=sign(y_tr).*acos(x_tr./sqrt(x_tr.^2+y_tr.^2))+2*pi;   % the azimuth angle

rho_x=lambda_min/2;
rho_y=c/2/B;

Nx=Na;Ny=Nf;

%%%% x_grid 和 y_grid 是设定的目标区域网格 %%%%%%%%%%%%%%%%%%%%%
%%%% x_range和y_range 是根据FT变换规则在N-点FT变换下的定标结果%%%%%%%%%%%%%%
%%%% 两者的关系是 y_grid 和 y_range
%%%% 本质是等价的，c/(f_step*Nf)=2rho_y,是不是就是说FT在距离维若采用fft变换只能分辨2倍rho_y的散射点?
%%%% 不是，因为有双程差因子在，factor*c/(f_step*Nf)=rho_y

%%%% x_range和x_grid的关系
%%%% 
factor=0.5;         % 双程差=0.5；单程差=1；
x_grid=[(-Nx/2:Nx/2-1).']*rho_x/2;
y_grid=[(-Ny/2:Ny/2-1).']*rho_y/2;
y_range=factor*((-Nf/2:Nf/2-1).')*c/f_step/Nf;  
x_range=factor*(((0:Na-1)-Na/2).')./x_step*abs(sqrt(y_tr.^2+x_tr.^2)')*c/max(f)*ones(Na,1)./Na./Na; 



OBJECT=[0,2,1
        6,0,1
        9,4,1 ]*diag([rho_x rho_y 1]);      % [x,y,amp]
    
% OBJECT=[0,0,1]*diag([rho_x rho_y 1]);      % [x,y,amp]  
    
Dx=max(OBJECT(:,1))-min(OBJECT(:,1));   % range of x direction
Dy=max(OBJECT(:,2))-min(OBJECT(:,2));   % range of y direction

theta_max=2*asin((Lx+Dx)/sqrt((2*R0-Dy)^2+(Lx+Dx)^2));  % beamwidth in radius


% figure,title('目标坐标系下的目标和TR位置')
% plot(x_tr,y_tr,'ro','MarkerSize',5);hold on, plot(OBJECT(:,1),OBJECT(:,2),'*','MarkerSize',5)
% grid on; 

ECHO=zeros(Na,Nf);     % the matrix of superposed echo
R_object=[];
j=sqrt(-1);
for array_index=1:Na
    s=zeros(1,Nf);      % the echo towards each sampling position 
    for object_index=1:size(OBJECT,1)
        object_x=OBJECT(object_index,1);
        object_y=OBJECT(object_index,2);
        object_amp=OBJECT(object_index,3);
        R=sqrt((x_tr(array_index)-object_x)^2+(object_y-y_tr(array_index))^2);
        R_object=[R_object;R];
        s=s+object_amp*exp(-j*2*pi*f*2*R/c);
    end
    ECHO(array_index,:)=s;     
end

figure,plot(f,real(s),f,imag(s),'r'),legend('real','imag'),title('received signal')
figure,imagesc(real(ECHO))
figure,surf(real(ECHO));

%% 2D imaging
S_xFT_MP=zeros(Na,Nf);
kr=2*pi*f/c;    % uniform sampling in kr
kr=kr(:);       % being a column

Kx=zeros(Na,Nf);
Ky=zeros(Na,Nf);
for f_index=1:Nf
    kx=2*kr(f_index)*x_tr./sqrt(x_tr.^2+y_tr.^2);
    ky=sqrt(4*kr(f_index).^2-kx.^2);
    Kx(:,f_index)=kx;
    Ky(:,f_index)=ky;
end
% figure,plot(Kx,Ky,'.'),axis([min(min(Kx)) max(max(Kx)) min(min(Ky)) max(max(Ky))]),grid on;
% title(['kxky 分布 kr范围[' num2str(min(2*kr)) ' ' num2str(max(2*kr)) ']']),xlabel('kx'),ylabel('ky');

%% 先对x维做傅里叶变换，e^(-jk_x*x)
for f_index=1:Nf 
    temp=ECHO(:,f_index);
    kx=Kx(:,f_index);
%     x_grid=[((0:Na-1)-Na/2).']./x_step.*abs(sqrt(y_tr.^2+x_tr.^2))*c/f(f_index)./Na;
    x_fft=exp(-j*x_tr*(kx'))*temp;   %%% the key is doing FFT on x-position of receiver
    S_xFT_MP(:,f_index)=x_fft;
end
% S_xFT_MP=fftshift(fft(ECHO));
PHASE_COMPEN=exp(j.*Ky.*R0);        % 影响偏移量
S_xFT_MP=S_xFT_MP.*PHASE_COMPEN;
% 
% figure,plot([real(x_fft) imag(x_fft)]);
% surf(abs(S_xFT_MP)); title('方位向FT+匹配滤波')
% xlabel('kx');ylabel('kr');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%% 方法一 忽略ky维的非均匀性，直接利用公式计算IFT，e^(jk_x*x+jk_y*y)  %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 证明直接用公式推导是正确合理的，这个分辨率较差，就是匹配滤波 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_yIFT=zeros(Na,Nf);
for array_index=1:Na
    temp=S_xFT_MP(array_index,:);
    ky=Ky(array_index,:).';
    y_ift=temp*exp(j.*ky*(y_grid.'))./Nf;
    S_yIFT(array_index,:)=y_ift;
end
% figure,surf(abs(S_yIFT)),title('1D sum');

S_xyIFT=zeros(Na,Nf);
for f_index=1:Nf
    temp=S_yIFT(:,f_index);
    kx=Kx(:,f_index);
    x_ift=exp(j.*x_grid*(kx.'))*temp./Na;
    S_xyIFT(:,f_index)=x_ift;
end
% figure
% surf(abs(S_xyIFT));

x_range=x_grid;y_range=y_grid;
[X,Y] = meshgrid(x_range,y_range);
figure,
% subplot(1,2,1),
surf(X,Y,abs(S_xyIFT).'),title('2D Sum');
figure
% subplot(1,2,2),
imagesc(x_range,y_range,mat2gray(abs((S_xyIFT).')));
xlabel('x(m)'),ylabel('y(m)');
hold on,plot(OBJECT(:,1),OBJECT(:,2),'wo','MarkerSize',6,'LineWidth',1.5);