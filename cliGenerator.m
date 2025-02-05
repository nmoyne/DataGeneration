% This script will generate 256 datapoints
% A rectangle of scatterers with a random M (in [1, 5]e-3)
% Between 1 and 3 circles with random centers, radius, and M perturbation between
% [1, 5]e-4 and 
% Final image saved as jpeg and terrain truth as json

l_pict = 15e-2;
R = 128;

% Ensure the output directory exists
outputDir = 'img';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Rectangle with x in [-7.5e-2, 7.5e-2] and z in [0, 15e-2]
rectangle = @(x, z) (x >= -l_pict/2 & x <= l_pict/2) & (z >= 0 & z <= l_pict);

% Pre-allocate a cell array to store results
temp_data = cell(16, 16);

parfor k = 1:64
    for j = 1:16 
        M1 = rand() * 4e-3 + 1e-3;
        n = randi([1, 3]);

        circles = struct([]);
        func_list = {rectangle};
        m_list = [M1];

        for i = 1:n
            x_c = (rand() - 0.5) * l_pict/2;
            z_c = rand() * l_pict/2;
            r_c = rand() * 3e-2;
            m_c = rand() * 9e-3 + 1e-3;

            circle = @(x, z) ((x - x_c)^2 + (z - z_c)^2) <= r_c^2;
            func_list{end + 1} = circle;
            m_list = [m_list, m_c];

            circles(i).x = x_c;
            circles(i).z = z_c;
            circles(i).r = r_c;
            circles(i).m = m_c;
        end

        rectangleEmpty = @(x, z) rectangle(x, z);
        for i = 2:length(func_list) % Skip the first rectangle function
            rectangleEmpty = @(x, z) rectangleEmpty(x, z) & ~func_list{i}(x, z);
        end

        I = uint8(ones(256, 256) * R);
        [xs, zs, RC] = fieldScatterers(func_list, m_list, l_pict, I);
        
        % Generate the ultrasound image
        [xi, zi, I_res] = buildPicture(l_pict, xs, zs, RC);

        % Save the image
        img_name = sprintf('%d_%d_bmode.jpeg', k, j);
        img_path = fullfile(outputDir, img_name);
        imwrite(I_res, img_path);

        % Store results in a cell array (instead of modifying terrain_data directly)
        temp_data{k, j} = struct('name', img_name, 'circles', circles, 'level', M1);

        % Print iteration completion message (unordered)
        fprintf('Iteration %d completed for parallel iteration %d\n', j, k);
    end
end

% Convert cell array to a linear struct array
terrain_data = [temp_data{:}];

% Save JSON file
jsonStr = jsonencode(terrain_data, 'PrettyPrint', true);
fid = fopen('terrain_data.json', 'w');
fwrite(fid, jsonStr, 'char');
fclose(fid);

disp('Processing complete. Images saved in "img/" and JSON saved as "terrain_data.json".');







