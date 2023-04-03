%% 仿真信号输入
clc
clear
close all

Fs = 20000;              % sampling frequency (in Hz)
T = 10;                  % time length (in s)
L = T*Fs;                % signal length (number of samples)
f0 = 0.01*Fs;            % cycle frequency (in Hz)
t = 0:1/Fs:(T-1/Fs);

%%%调制信号设定%%%
Am1=1;
fm1=5;
Am2=1;
fm2=10;
Am3=1;
fm3=15;
Am4=1;
fm4=20;
ym1=Am1*cos(2*pi*fm1*t);
ym2=Am2*cos(2*pi*fm2*t);
ym3=Am3*cos(2*pi*fm3*t);
ym4=Am4*cos(2*pi*fm4*t);
%%%宽带载波设定%%%
y_white=wgn(1,L,40000,'linear');
fs=Fs/2;
Wp=[2500 7500]/fs;
Ws=[2000 8000]/fs;
Rp=3;
Rs=40;
[n,Wn]=buttord(Wp,Ws,Rp,Rs);
[b,a]=butter(n,Wn);
yc=filter(b,a,y_white);
%%生成总体信号%%%
y=(1+ym1+ym2+ym3+ym4).*yc;
snr=-20:0;

pjjt=zeros(1,4);
pd=zeros(length(snr),4);
xh=1000;%循环次数
cunchu=zeros(length(snr),4,xh);
for is=1:length(snr)
for sy=1:xh
x=awgn(y,snr(is),'measured');
% x=awgn(y,-5,'measured');
% y=x;

% figure
% plot(t,y)
% title('Synthetic cyclostationary signal')
% xlabel('time (s)')
Nw = 64;               % window length (number of samples),df ~ 1.5*Fs/Nw (in Hz)
alpha_max =50;       % maximum cyclic frequency to scan (in Hz)
%% 快速循环平稳分析
                        % (should cover the cyclic frequency range of interest)
opt.coh = 1;            % compute sepctral coherence? (yes=1, no=0)

[S,alpha,f,Nv] = Fast_SC(x,Nw,alpha_max,Fs,opt);
%% 加权函数处理
SS=zeros(size(S));%自相关切片
w=zeros(size(S,1),1);%自相关加权函数
for i=1:length(f)
slice=abs(S(i,1:end));
sliceC=xcorr(slice);
rn=(length(sliceC)+1)/2;
sliceC=sliceC(rn:end);
SS(i,:)=sliceC;
rpeak=findpeaks(sliceC(2:end));
rmax=max(rpeak);
w(i)=(rmax/(sliceC(1)-rmax)).^2;
end

k(1,:)=w;%切片自相关加权函数

YY = fft(S');
k(2,:) = kurtosis(abs(YY));%切片傅里叶峭度

k=k./repmat(max(k,[],2),1,length(k));%归一化
Sk=zeros(size(k,1),length(alpha)-1);%加权EES
for i=1:size(k,1)
    kc=repmat(k(i,:)',1,size(S,2)-1);
    if i==2
        kcs=mean(abs((S(:,2:end).*kc).^2));
    else
        kcs=mean(abs(S(:,2:end).*kc));
    end
    Sk(i,:)=kcs;
end
Sk(i+1,:)=mean(abs(S(:,2:end)));%EES

nlevel = 4;     % number of decomposition levels
prewhiten = 1;
if prewhiten == 1
   x = x - mean(x);
   Na = 100;
   a = lpc(x,Na);
   x = fftfilt(a,x);
   x = x(Na+1:end);		% it is very important to remove the transient of the whitening filter, otherwise the SK will detect it!!
end
[~,~,~,~,~,~,bw,fc]= Fast_kurtogram2(x,nlevel,Fs);
[~,f1]=min(abs(f-(fc-bw/2)));
[~,f2]=min(abs(f-(fc+bw/2)));
Sk(i+2,:)=mean(abs(S(f1:f2,2:end)));%快速峭度EES
%% 性能评价指标
%主要经验参数
Nthw =40;%平滑窗长
Nthv = fix(0.9*Nthw);%重叠率
level=0.9;%排序等级
X=0.01;%带宽系数

%排序法阈值曲线
Nth=ceil(length(Sk)/(Nthw-Nthv));%窗口数量
thi=zeros(size(Sk,1),Nth);%窗口阈值曲线
ath=zeros(1,Nth);%窗口对应循环频率
for i=1:Nth
    N1=(Nthw-Nthv)*(i-1)+1;
    N2=min(((Nthw-Nthv)*(i-1)+Nthw),length(Sk));
    thic=sort(Sk(:,N1:N2),2);%升序重排列窗口
    thin=round(level*(N2-N1+1));%窗口阈值索引
    thi(:,i)=thic(:,thin);
    ath(i)=alpha(ceil((N1+N2)/2));
end
th= interp1(ath,thi',alpha(2:end),'pchip','extrap');%插值生成阈值曲线
th=1*th';%添加系数

%检测带宽
ac=[fm1;fm2;fm3;fm4;25;30;35;40];%各个检测频带中心频率
B=zeros(2*length(ac),1);%检测带宽
Bn=zeros(2*length(ac),1);%检测带宽序号
m=zeros(size(Sk,1),length(ac));%检测带宽内峰值
mn=zeros(size(Sk,1),length(ac));%检测带宽内峰值序号
for bi=1:length(ac)
    B(2*bi-1)=ac(bi)-X*ac(1);
    B(2*bi)=ac(bi)+X*ac(1);
    [~,Bn(2*bi-1)]=min(abs(alpha(2:end)-B(2*bi-1)));
    [~,Bn(2*bi)]=min(abs(alpha(2:end)-B(2*bi)));
    [m(:,bi),mn(:,bi)]=max(Sk(:,Bn(2*bi-1):Bn(2*bi)),[],2);
    mn(:,bi)=mn(:,bi)+Bn(2*bi-1)-1;
end
col = mn';
col=col(:);
row=zeros(length(col),1);
for i=1:size(Sk,1)
    row(1+(i-1)*length(ac):i*length(ac)) = i.*ones(length(ac),1);
end
ind = sub2ind(size(th),row,col);
thn=reshape(th(ind),length(ac),size(Sk,1))';%检测带宽峰值对应阈值

pm=(m./thn)';%峰值和阈值比
ppm=geomean(pm);%几何平均峰值比
% ppm=mean(pm);
% ppm=harmmean(pm);
cunchu(is,:,sy)=ppm;%存储数组
pjjt=pjjt+ppm;
end
pd(is,:)=pjjt/xh;
pjjt=zeros(1,4);
end
%% 画图
color1=[209,42,32;21,29,41;80,138,178;250,192,61];
color2=color1/255;
saveIndex=0;
%ATR值
ATR=mean(cunchu,3);
exl=[snr' ATR];
exl=exl(1:end,:);
figure
set(gcf,'unit','centimeters','position',[10 5 18 8]); % 10cm*17.4cm
%     set(gcf,'ToolBar','none','ReSize','off');   % 移除工具栏
set(gcf,'color','w'); % 背景设为白色
colororder(color2)
plot(exl(:,1),exl(:,2),'s-','LineWidth',1.2)
hold on 
plot(exl(:,1),exl(:,3),'o--','LineWidth',1.2)
plot(exl(:,1),exl(:,4),'x:','LineWidth',1.2)
plot(exl(:,1),exl(:,5),'d-.','LineWidth',1.2)
grid on
xlabel('SNR[dB]'),ylabel('ATR')
legend('WEES','AWES','EES','Kurtogram')
legend(Location="northwest",Box="off");
if saveIndex==1
    exportgraphics(gca, 'ATR.tiff')
end

%单个准确率图
npdd=mean(cunchu>1.5,3);
snr=snr(1:end);
npdd=npdd(1:end,:);
figure
set(gcf,'unit','centimeters','position',[10 5 18 8]); % 10cm*17.4cm
%     set(gcf,'ToolBar','none','ReSize','off');   % 移除工具栏
set(gcf,'color','w'); % 背景设为白色
set(gca,"ColorOrder",color2);
plot(snr,npdd(:,1),'s-',snr,npdd(:,2),'o--',snr,npdd(:,3),'x:',snr,npdd(:,4),'d-.','LineWidth',1.2)
colororder(color2)
ylim([0,1.1])
grid on
xlabel('SNR[dB]'),ylabel('Detection probability')
legend('WEES','AWES','EES','Kurtogram')
legend(Location="southeast",Box="off");
if saveIndex==1
    exportgraphics(gca, 'T=1.5.tiff')
end

npdd=mean(cunchu>2,3);
snr=snr(1:end);
npdd=npdd(1:end,:);
figure
set(gcf,'unit','centimeters','position',[10 5 18 8]); % 10cm*17.4cm
%     set(gcf,'ToolBar','none','ReSize','off');   % 移除工具栏
set(gcf,'color','w'); % 背景设为白色
set(gca,"ColorOrder",color2);
plot(snr,npdd(:,1),'s-',snr,npdd(:,2),'o--',snr,npdd(:,3),'x:',snr,npdd(:,4),'d-.','LineWidth',1.2)
colororder(color2)
ylim([0,1.1])
grid on
xlabel('SNR[dB]'),ylabel('Detection probability')
legend('WEES','AWES','EES','Kurtogram')
legend(Location="southeast",Box="off");
if saveIndex==1
    exportgraphics(gca, 'T=2.tiff')
end

%不同T下的准确率
cunzaiyz=1:0.1:3;
bijiao=cunchu(6,:,:);
zql=zeros(length(cunzaiyz),4);
for i=1:length(cunzaiyz)
zql(i,:)=mean(bijiao>cunzaiyz(i),3);
end
figure
set(gcf,'unit','centimeters','position',[10 5 18 8]); % 10cm*17.4cm
%     set(gcf,'ToolBar','none','ReSize','off');   % 移除工具栏
set(gcf,'color','w'); % 背景设为白色
set(gca,"ColorOrder",color2);
plot(cunzaiyz,zql(:,1),'s-',cunzaiyz,zql(:,2),'o--',cunzaiyz,zql(:,3),'x:',cunzaiyz,zql(:,4),'d-.','LineWidth',1.2)
colororder(color2)
ylim([0,1.1])
grid on
xlabel('Threshold'),ylabel('Detection probability')
legend('WEES','AWES','EES','Kurtogram')
legend(Box="off");
if saveIndex==1
    exportgraphics(gca, '-15dB.tiff')
end