function [xs, zs, RC] = specialScatterers(type, R, l_pict)
    switch type
        case 'one'
            % One scatterer at the center
            xs = 0;
            zs = (l_pict / 4);
            RC = R;
        
        case 'arc'
            % 5 scatterers on an arc from -30° to 30° at distance l_pict/2
            angles = linspace(-30, 30, 5) * pi / 180; % Convert to radians
            radius = l_pict / 2;
            xs = radius * sin(angles);
            zs = radius * cos(angles);
            RC = R * ones(size(xs)); % Same reflection coefficient for all
        
        case 'arc_mult'
            % 3 sets of 5 scatterers on arcs at different distances
            distances = [l_pict / 4, l_pict / 3, l_pict / 2]; % Arc radii
            angles = linspace(-30, 30, 5) * pi / 180; % Convert to radians
            
            xs = [];
            zs = [];
            RC = [];
            
            for radius = distances
                xs = [xs, radius * sin(angles)];
                zs = [zs, radius * cos(angles)];
                RC = [RC, R * ones(size(angles))]; % Assign same RC for each set
            end
        
        case 'line'
            zs = [2, 3, 5, 7, 10]*1e-2;
            xs = zeros(size(zs));
            RC = R*ones(size(zs));

            
        otherwise
            error('Invalid type. Choose "one", "arc", "line" or "arc_mult".');
    end
end
