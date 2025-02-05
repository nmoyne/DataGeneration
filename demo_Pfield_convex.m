clc;clear;close all;

param = getparam('C5-2v');

R = param.radius;
p = param.pitch;
N = param.Nelements;
L = 2*R*sin(asin(p/2/R)*(N-1)); % chord length
h = sqrt(R^2-L^2/4);
th = linspace(atan2(-L/2,h),atan2(L/2,h),N);
xe = R*sin(th);
ze = R*cos(th)-h;

param_suba = param;
param_suba.Nelements = 32;

xf = [0]; zf = 30e-3; % focus position (in m)

txdel_suba = txdelay(xf,zf,param_suba); % in s

xs=[-20;0;0;0;0;0]/1000;zs=[10;20;30;40;50;75]/1000;RC=ones(numel(zs),1);
param.fs=4*param.fc;

I = imread('heart.jpg');
I = rgb2gray(I);
[xs,~,zs,RC] = genscat([NaN 15e-2],param,I);

option.WaitBar = false;
for i=1:97
    subaper=(1:32)+(i-1);
    txdel_1to32 = NaN(1,128);
    txdel_1to32(subaper) = txdel_suba;
    RF{i}=simus(xs,zs,RC,txdel_1to32,param,option);
    txdel{i}=txdel_1to32;
end

% IQ demodulation
IQ = cell(numel(RF),1);  % this cell will contain the I/Q series

for k = 1:numel(IQ)
    IQ{k} = rf2iq(RF{k},param.fs,param.fc);
end

% Formation des voies 
param.fnumber = [];
[xi,zi] = impolgrid([256 256],90e-3,param);

bIQ = zeros(256,256,numel(RF));  % this array will contain the 13 images

h = waitbar(0,'Démodulation DAS');
for k = 1:numel(IQ)
    waitbar(k/numel(IQ),h,['DAS: I/Q series #' int2str(k) ' of 13'])
    bIQ(:,:,k) = das(IQ{k},xi,zi,txdel{k},param);
end
close(h)

% Time Gain compensation
%bIQ = tgc(bIQ);

% Formation de l'image en BMODE par somme des images
bIQ_sum = sum(bIQ, 3);
I = bmode(bIQ_sum(:,:,1),40); % log-compressed image

% Display the final B-mode image
figure;
pcolor(xi*100,zi*100,I)
c = colorbar;
c.YTick = [0 255];
c.YTickLabel = {'-log dB','0 dB'};
colormap gray
title('Simulated ultrasound image')
ylabel('[cm]')
shading interp
axis equal ij tight
set(gca,'XColor','none','box','off')