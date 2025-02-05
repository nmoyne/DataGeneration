% Scatterers : must be drawn from two distributions.

% Uniform between 0 and 4e-2 (x, z)
% Plus a line on the diagonal

xs1 = rand(1, 5000)*8*1e-2 - 4*1e-2; % in m
zs1 = rand(1, 5000)*4*1e-2;

xs2 = rand(1, 1000)*8*1e-2;
ds2 = rand(1, 1000)*1e-3;
zs2 = xs2 + ds2;

xs = [xs1, xs2];
zs = [zs1, zs2];

RC = ones(size(xs));

% TO BE CHANGED
param = getparam('L11-5v');

% Depends on params ?
tilt = linspace(-10,10,21)/180*pi; % tilt angles in rad
txdel = cell(21,1); % this cell will contain the transmit delays

for k = 1:21
    txdel{k} = txdelay(param,tilt(k));
end

RF = cell(17,1); % this cell will contain the RF series
param.fs = 4*param.fc; % sampling frequency in Hz

option.WaitBar = false; % remove the progress bar of SIMUS
for k = 1:21
    RF{k} = simus(xs,zs,RC,txdel{k},param,option);
end

IQ = cell(21,1);  % this cell will contain the I/Q series

for k = 1:21
    IQ{k} = rf2iq(RF{k},param.fs,param.fc);
end

param.fnumber = [];

[xi,zi] = meshgrid(linspace(-2e-2,2e-2,200),linspace(0,4e-2,200));

bIQ = zeros(200,200,21);  % this array will contain the 21 I/Q images

h = waitbar(0,'');
for k = 1:21
    waitbar(k/21,h,['DAS: I/Q series #' int2str(k) ' of 21'])
    bIQ(:,:,k) = das(IQ{k},xi,zi,txdel{k},param);
end
close(h)

cIQ = sum(bIQ,3); % this is the compound beamformed I/Q
I = bmode(cIQ,40); % log-compressed image
imagesc(xi(1,:)*1e2,zi(:,1)*1e2,I)
colormap gray
title('Compound PW-based echo image')

axis equal ij
set(gca,'XColor','none','box','off')
c = colorbar;
c.YTick = [0 255];
c.YTickLabel = {'-40 dB','0 dB'};
ylabel('[cm]')
