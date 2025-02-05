function [xs, zs, RC] = generateScatterers(R, M1, M2, l_pict, l_scat, x_lim, z_lim)
    % Define image size
    ImgSize = 256;

    % Create the uniform gray image
    I = uint8(ones(ImgSize, ImgSize) * R);

    % Create 2D distribution of scatterers
    [xs1, ~, zs1, RC1] = genscat([NaN, l_pict], M1, I);
    [xs2, ~, zs2, RC2] = genscat([NaN, l_scat], M2, I);

    % Find indices of scatterers in the first set that fall within the second grid
    idx_to_remove = (xs1 >= x_lim & xs1 <= x_lim + l_scat) & ...
                    (zs1 >= z_lim & zs1 <= z_lim + l_scat);

    % Construct the final lists:
    % - All scatterers from the second grid are included
    % - All scatterers from the first grid are included except those in the overlap region
    xs = [xs1(~idx_to_remove); xs2 + x_lim + l_scat/2];
    zs = [zs1(~idx_to_remove); zs2 + z_lim];
    RC = [RC1(~idx_to_remove); RC2];
end


















