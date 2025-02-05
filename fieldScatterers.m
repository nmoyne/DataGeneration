function [xs, zs, RC] = fieldScatterers(distr_list, m_list, l_pict, I)

    xs = [];
    zs = [];
    RC = [];

    for k = 1:numel(distr_list)
        [xsk, ~, zsk, RCk] = genscat([NaN, l_pict], m_list(k), I);
        filter = distr_list{k};

        mask = arrayfun(@(x, z) filter(x, z), xsk, zsk);
        
        % Append only valid values
        xs = [xs; xsk(mask)];
        zs = [zs; zsk(mask)];
        RC = [RC; RCk(mask)];
    end
end