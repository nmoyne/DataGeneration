% Define parameters
R = 128;
M1 = 4e-3;
M2 = 1e-3;
l_pict = 15e-2;
l_scat = 2e-2;
x_lim = -2e-2;
z_lim = 5e-2;

% Generate scatterers
%[xs, zs, RC] = generateScatterers(R, M1, M2, l_pict, l_scat, x_lim, z_lim);
[xs, zs, RC] = specialScatterers("one", R, l_pict);

% Circle centered at (0, 4e-2) with radius 1
circle1 = @(x, z) (x.^2 + (z - 4e-2).^2) <= (1e-2)^2;

% Circle centered at (3e-2, 6e-2) with radius 1
circle2 = @(x, z) ((x - 3e-2).^2 + (z - 6e-2).^2) <= (1e-2)^2;

% Circle centered at (-3e-2, 6e-2) with radius 1
circle3 = @(x, z) ((x + 3e-2).^2 + (z - 6e-2).^2) <= (1e-2)^2;

% Rectangle with x in [-7.5e-2, 7.5e-2] and z in [0, 15e-2]
rectangle = @(x, z) (x >= -7.5e-2 & x <= 7.5e-2) & (z >= 0 & z <= 15e-2);

% Rectangle without the circles
rectangleEmpty = @(x, y) rectangle(x, y) ...
& ~circle3(x,y) ...
& ~circle2(x,y) ...
& ~circle1(x,y);


func_list = {rectangleEmpty, circle2, circle3};
m_list = [4, 1, 1]*1e-3;
I = uint8(ones(256, 256) * R);

%[xs, zs, RC] = fieldScatterers(func_list, m_list, l_pict, I);


figure;
scatter(xs*1e2,zs*1e2,20,RC,'filled')
colormap([hot; 1-hot])
axis equal ij tight
set(gca,'XColor','none','box','off')
title([int2str(numel(RC)) ' scatterers'])
ylabel('[cm]')

% Call the function to generate the ultrasound image
[xi, zi, I] = focusedBuildPicture(l_pict, xs, zs, RC);

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