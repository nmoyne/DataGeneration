% Request image name from the user
img_name = input('Enter the image name (e.g., "1_1_bmode.jpeg"): ', 's');

% Load JSON file
json_file = 'terrain_data.json';
if exist(json_file, 'file') ~= 2
    error('JSON file not found: %s', json_file);
end

fid = fopen(json_file, 'r');
raw = fread(fid, inf, 'uint8=>char')';
fclose(fid);
terrain_data = jsondecode(raw);

% Find the corresponding entry in JSON
img_data = [];
for i = 1:length(terrain_data)
    if strcmp(terrain_data(i).name, img_name)
        img_data = terrain_data(i);
        break;
    end
end

% Check if image data was found
if isempty(img_data)
    error('Image not found in terrain_data.json');
end

% Load and display image
img_path = fullfile('img', img_name);
if exist(img_path, 'file') ~= 2
    error('Image file not found: %s', img_path);
end

I = imread(img_path);
figure;
imshow(I);
hold on;

% Title: "Base MAINDIST : level"
title(sprintf('Base MAINDIST : %.4f', img_data.level), 'FontSize', 14, 'FontWeight', 'bold');

% Draw circles and display perturbation
num_circles = length(img_data.circles);
legendEntries = {};
for i = 1:num_circles
    % Get circle parameters
    x = img_data.circles(i).x;
    z = img_data.circles(i).z;
    r = img_data.circles(i).r;
    m = img_data.circles(i).m;
    
    % Compute perturbation
    perturbation = m - img_data.level;
    
    % Convert coordinates for image display (assuming standard Cartesian)
    x_pixel = round((x + 7.5e-2) / (15e-2) * size(I, 2));
    z_pixel = round(z / (15e-2) * size(I, 1));

    % Draw red circle
    theta = linspace(0, 2*pi, 100);
    xc = x_pixel + r * cos(theta) * size(I, 2) / (15e-2);
    zc = z_pixel + r * sin(theta) * size(I, 1) / (15e-2);
    plot(xc, zc, 'r', 'LineWidth', 2);

    % Display perturbation inside the circle if space allows
    if num_circles <= 3
        text(x_pixel, z_pixel, sprintf('%.2e', perturbation), ...
            'Color', 'red', 'FontSize', 10, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center');
    else
        % Store legend entry if too many circles
        legendEntries{end+1} = sprintf('Perturbation %.2e at (%.3f, %.3f)', perturbation, x, z);
    end
end

% Display legend if too many circles
if num_circles > 3
    legend(legendEntries, 'TextColor', 'red', 'FontSize', 10, 'Location', 'southoutside');
end

hold off;
disp('Image displayed with annotations.');
