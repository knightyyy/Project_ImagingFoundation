%% ��άƽ������㷨�����ڸ���Ҷ�任
clc;
clear all;
close all;

%% ��ʼ�����趨
c=3e8;%����
f0=34e9;%��ʼƵ��
B=20e9;%����
lambda_min=c/(f0+B/2);%����


%% ��ά�ֱ���ȷ��
XResolution=lambda_min/4;%����ά�ֱ���
YResolution=lambda_min/4;%�߶�ά�ֱ���
ZResolution=c/(2*B);%����ά�ֱ���

%% �۲��������
Nx=250;%x����Ĳ�������
SampleIntervalX=lambda_min/4;%x����Ĳ������
PositionDetX=(-Nx/2:Nx/2-1)*SampleIntervalX;%x����������λ��
Ny=250;%y����Ĳ�������
SampleIntervalY=lambda_min/4;%y����Ĳ������
PositionDetY=(-Ny/2:Ny/2-1)*SampleIntervalY;%y����������λ��



%% ����ά���
Nf=100;%����ά��������
rangeZ=ZResolution*Nf;
% SampleIntervalf=B/Nf;%����ά�������
SampleIntervalf=c/(2*rangeZ);%����ά�������
f=f0+(0:Nf-1)*SampleIntervalf;%����Ƶ��
k=2*pi*f/c;%������Ƶ��
PositionDetZ=0;%���ù۲�ƽ��ľ���άλ��
%% Ŀ��ƽ�����ɲⷶΧ;�������ߵķ���ǽ�С�����ƵĿ�����Ϊ��һ����
rangeX=max(PositionDetX)-min(PositionDetX);
rangeY=max(PositionDetY)-min(PositionDetY);

%% Ŀ��ƽ����������
N_xgrid=Nx;
N_ygrid=Ny;%������Ŀ�����������۲�ƽ����ͬ
N_zgrid=Nf;
IntervalXgrid=rangeX/N_xgrid;
IntervalYgrid=rangeY/N_ygrid;
IntervalZgrid=rangeZ/N_zgrid;
x_grid=(-N_xgrid/2:N_xgrid/2-1)*IntervalXgrid;
y_grid=(-N_ygrid/2:N_ygrid/2-1)*IntervalYgrid;
z_grid=(-N_zgrid/2:N_zgrid/2-1)*IntervalZgrid;

%% Ŀ���λ������
% val=diag([XResolution,YResolution,ZResolution]);
Ptar=[10 10 5
    -30 -30 20
50 50 5
]*diag([XResolution,YResolution,ZResolution]);
Ntar=length(Ptar(:,1));%����Ŀ��ĸ���

%% �ز��źŹ���
S_echo=zeros(Nx,Ny,Nf);%�����ûز��źž���
for index_x=1:Nx;
    for index_y=1:Ny;
            s=0;
            for index_tar=1:Ntar;
                xtar=Ptar(index_tar,1);
                ytar=Ptar(index_tar,2);
                ztar=Ptar(index_tar,3);
                sigma=1;%��Ŀ���ķ�ֵΪ1��
                R=sqrt((xtar-PositionDetX(index_x)).^2+(ytar-PositionDetY(index_y)).^2+(ztar-PositionDetZ).^2);
                s=s+sigma*exp(-1j*2.*k.*R);
            end
            S_echo(index_x,index_y,:)=s;
    end
end

%% ������Ƶ�ʼ���
% for index_f=1:Nf;
%     for index_y=1:Ny;
%         Kx_range=2*pi*f(index_f)*rangeX/(c*sqrt((rangeX/2).^2+(rangeZ).^2));
%         KX1=linspace(-Kx_range,Kx_range,Nx);
%         Kx(:,index_y,index_f)=KX1;
%     end
%     for index_x=1:Nx;
%         Ky_range=2*pi*f(index_f)*rangeY/(c*sqrt((rangeY/2).^2+(rangeZ).^2));
%         KY1=linspace(-Ky_range,Ky_range,Ny);
%         Ky(index_x,:,index_f)=KY1;
%     end
% end
for index_f=1:Nf;
    for index_y=1:Ny;
        Kx_range=2*pi*f(index_f)*(PositionDetX)./(c*sqrt((PositionDetX).^2+(rangeZ).^2));
        Kx(:,index_y,index_f)=Kx_range;
    end
    for index_x=1:Nx;
        Ky_range=2*pi*f(index_f)*(PositionDetY)./(c*sqrt((PositionDetY).^2+(rangeZ).^2));
        Ky(index_x,:,index_f)=Ky_range;
    end
end
for index_x=1:Nx;
    for index_y=1:Ny;
%         Kx1=zeros(1,Nf);
%         Ky1=zeros(1,Nf);
%         Kx1=reshape(Kx(index_x,index_y,:),1,Nf);
%         Ky1=reshape(Ky(index_x,index_y,:),1,Nf);
%         Kz(index_x,index_y,:)=sqrt((4*k.^2)-Kx1.^2-Ky1.^2);
        Kz(index_x,index_y,:)=2*k;
    end
end

%% ��ά����Ҷ�任�����õķ�ʽ�ǹ̶�Ƶ��ά��֮��Զ�ά��λά���и���Ҷ�任
for index_f=1:Nf;
    KX(:,:)=Kx(:,:,index_f);
    S_echoX(:,:)=S_echo(:,:,index_f);
    for index_y=1:Ny;
        S_FTx(:,index_y)=exp(-1j*PositionDetX.'*KX(:,index_y).')*S_echoX(:,index_y);
    end
    KY(:,:)=Ky(:,:,index_f);
    for index_x=1:Nx;
        S_FTy(index_x,:)=exp(-1j*PositionDetY.'*KY(index_x,:))*(S_FTx(index_x,:).');
    end
    S_FTxy1(:,:,index_f)=S_FTy;
end

%% ��λ����
for index_x=1:Nx;
    for index_y=1:Ny;
        R0=PositionDetZ;
        Kz1=Kz(index_x,index_y,:);
        Phase_comp=exp(1j.*Kz1.*R0);
        S_FTxy1(index_x,index_y,:)=S_FTxy1(index_x,index_y,:).*Phase_comp;
    end
end
%% kz������Ƶ�ʲ�ֵ���Ȼ�����
% for index_x=1:Nx;
%     KZ1(:,:)=Kz(index_x,:,:);
%     for  index_y=1:Ny;
%         KZ2=KZ1(index_y,:);
%         KZmax=max(KZ2);
%         Kz(index_x,index_y,:)=linspace(-KZmax,KZmax,Nf);
%     end
% end

%% ��ά����Ҷ��任
for index_x=1:Nx;
    KZ(:,:)=Kz(index_x,:,:);
    S_bIFTz(:,:)=S_FTxy1(index_x,:,:);
    for index_y=1:Ny;
        S_IFTz(index_y,:)=exp(1j*KZ(index_y,:).'*z_grid )*S_bIFTz(index_y,:).';
    end
    S_IFTxyz(index_x,:,:)=S_IFTz;
end
for index_f=1:Nf;
    KX(:,:)=Kx(:,:,index_f);
    S_comp(:,:)=S_IFTxyz(:,:,index_f);
    for index_y=1:Ny;
        S_IFTx(:,index_y)=exp(1j*KX(:,index_y)*x_grid)*S_comp(:,index_y);
    end
    KY(:,:)=Ky(:,:,index_f);
    for index_x=1:Nx;
        S_IFTy(index_x,:)=exp(1j*KY(index_x,:).'*y_grid)*(S_IFTx(index_x,:).');
    end
    S_IFTxy(:,:,index_f)=S_IFTy;
end

%% ���ø�Ƶ��ͼ���ںϵķ�ʽ����ͼ����֮��ͼ���ںϵķ�ʽ��Ϊfunction���Ӳ�ͬ���ò�ͬ�����ϵ���Ƭ��ʵ��Ŀ�����ά����
yxy=S_IFTxy;
[line,array]=size(yxy(:,:,1));
    for i=1:line;
        for j=1:array;
            y_f(i,j)=(sum(yxy(i,j,:)))./Nf;
        end
    end
yxz=S_IFTxy;
    for i=1:Nx;
        for j=1:Nf;
            y_xz(i,j)=(sum(yxz(i,:,j)))./Ny;
        end
    end
yyz=S_IFTxy;
for i=1:Ny;
    for j=1:Nf;
        y_yz(i,j)=(sum(yyz(:,i,j)))./Nx;
    end
end
% [LL,HL,LH,HH]=mydwt2(y1);
% y_f=myidwt2(LL,HL,LH,HH);

% S1=S_IFTxy;
% y_f=imagefusion(Nx,Ny,Nf,S1);
% y_xz=imagefusionxz(Nx,Ny,Nf,S1);

%% ��ͼ
figure
imagesc(x_grid,y_grid,mat2gray(abs((y_f).')))
hold on 
plot(Ptar(:,1),Ptar(:,2),'ro');
figure
imagesc(x_grid,z_grid,mat2gray(abs(y_xz).'))

hold on
plot(Ptar(:,1),Ptar(:,3),'ro');
% axis([min(x_grid) max(x_grid) min(x_grid) max(x_grid)]);
figure
imagesc(z_grid,y_grid,mat2gray(abs(y_yz)))
hold on
plot(Ptar(:,3),Ptar(:,2),'ro');
% axis([min(x_grid) max(x_grid) min(x_grid) max(x_grid)]);