%% �����ź�����
clc
clear
close all

Fs = 20000;              % sampling frequency (in Hz)
T = 10;                  % time length (in s)
L = T*Fs;                % signal length (number of samples)
f0 = 0.01*Fs;            % cycle frequency (in Hz)
t = 0:1/Fs:(T-1/Fs);

%%%�����ź��趨%%%
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
%%%����ز��趨%%%
y_white=wgn(1,L,40000,'linear');
fs=Fs/2;
Wp=[2500 7500]/fs;
Ws=[2000 8000]/fs;
Rp=3;
Rs=40;
[n,Wn]=buttord(Wp,Ws,Rp,Rs);
[b,a]=butter(n,Wn);
yc=filter(b,a,y_white);
%%���������ź�%%%
y=(1+ym1+ym2+ym3+ym4).*yc;
snr=-20:0;

pjjt=zeros(1,4);
pd=zeros(length(snr),4);
xh=1000;%ѭ������
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
%% ����ѭ��ƽ�ȷ���
                        % (should cover the cyclic frequency range of interest)
opt.coh = 1;            % compute sepctral coherence? (yes=1, no=0)

[S,alpha,f,Nv] = Fast_SC(x,Nw,alpha_max,Fs,opt);
%% ��Ȩ��������
SS=zeros(size(S));%�������Ƭ
w=zeros(size(S,1),1);%����ؼ�Ȩ����
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

k(1,:)=w;%��Ƭ����ؼ�Ȩ����

YY = fft(S');
k(2,:) = kurtosis(abs(YY));%��Ƭ����Ҷ�Ͷ�

k=k./repmat(max(k,[],2),1,length(k));%��һ��
Sk=zeros(size(k,1),length(alpha)-1);%��ȨEES
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
Sk(i+2,:)=mean(abs(S(f1:f2,2:end)));%�����Ͷ�EES
%% ��������ָ��
%��Ҫ�������
Nthw =40;%ƽ������
Nthv = fix(0.9*Nthw);%�ص���
level=0.9;%����ȼ�
X=0.01;%����ϵ��

%������ֵ����
Nth=ceil(length(Sk)/(Nthw-Nthv));%��������
thi=zeros(size(Sk,1),Nth);%������ֵ����
ath=zeros(1,Nth);%���ڶ�Ӧѭ��Ƶ��
for i=1:Nth
    N1=(Nthw-Nthv)*(i-1)+1;
    N2=min(((Nthw-Nthv)*(i-1)+Nthw),length(Sk));
    thic=sort(Sk(:,N1:N2),2);%���������д���
    thin=round(level*(N2-N1+1));%������ֵ����
    thi(:,i)=thic(:,thin);
    ath(i)=alpha(ceil((N1+N2)/2));
end
th= interp1(ath,thi',alpha(2:end),'pchip','extrap');%��ֵ������ֵ����
th=1*th';%���ϵ��

%������
ac=[fm1;fm2;fm3;fm4;25;30;35;40];%�������Ƶ������Ƶ��
B=zeros(2*length(ac),1);%������
Bn=zeros(2*length(ac),1);%���������
m=zeros(size(Sk,1),length(ac));%�������ڷ�ֵ
mn=zeros(size(Sk,1),length(ac));%�������ڷ�ֵ���
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
thn=reshape(th(ind),length(ac),size(Sk,1))';%�������ֵ��Ӧ��ֵ

pm=(m./thn)';%��ֵ����ֵ��
ppm=geomean(pm);%����ƽ����ֵ��
% ppm=mean(pm);
% ppm=harmmean(pm);
cunchu(is,:,sy)=ppm;%�洢����
pjjt=pjjt+ppm;
end
pd(is,:)=pjjt/xh;
pjjt=zeros(1,4);
end
%% ��ͼ
color1=[209,42,32;21,29,41;80,138,178;250,192,61];
color2=color1/255;
saveIndex=0;
%ATRֵ
ATR=mean(cunchu,3);
exl=[snr' ATR];
exl=exl(1:end,:);
figure
set(gcf,'unit','centimeters','position',[10 5 18 8]); % 10cm*17.4cm
%     set(gcf,'ToolBar','none','ReSize','off');   % �Ƴ�������
set(gcf,'color','w'); % ������Ϊ��ɫ
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

%����׼ȷ��ͼ
npdd=mean(cunchu>1.5,3);
snr=snr(1:end);
npdd=npdd(1:end,:);
figure
set(gcf,'unit','centimeters','position',[10 5 18 8]); % 10cm*17.4cm
%     set(gcf,'ToolBar','none','ReSize','off');   % �Ƴ�������
set(gcf,'color','w'); % ������Ϊ��ɫ
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
%     set(gcf,'ToolBar','none','ReSize','off');   % �Ƴ�������
set(gcf,'color','w'); % ������Ϊ��ɫ
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

%��ͬT�µ�׼ȷ��
cunzaiyz=1:0.1:3;
bijiao=cunchu(6,:,:);
zql=zeros(length(cunzaiyz),4);
for i=1:length(cunzaiyz)
zql(i,:)=mean(bijiao>cunzaiyz(i),3);
end
figure
set(gcf,'unit','centimeters','position',[10 5 18 8]); % 10cm*17.4cm
%     set(gcf,'ToolBar','none','ReSize','off');   % �Ƴ�������
set(gcf,'color','w'); % ������Ϊ��ɫ
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