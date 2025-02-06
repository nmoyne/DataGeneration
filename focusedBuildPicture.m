function [xi, zi, I] = buildPicture(l_pict, xs, zs, RC) 
    
    % A convex array
    param = getparam("C5-2v");
    
    % Calcul des délais pour un sub-array de 32
    param_suba = param;
    param_suba.Nelements = 32;

    xf = 0;
    zf = 7.5e-2;

    txdel_suba = txdelay(xf, zf, param_suba);

    % Simulation du signal RF

    RF = cell(13,1); % Stackage des RF
    txdel = cell(numel(RF), 1); % Délais (utiles pour beamforming)
    param.fs = 4*param.fc; % sampling frequency in Hz

    option.WaitBar = false;
    %h = waitbar(0,'Génération RF');
    for k = 1:numel(RF)
        %waitbar(k/numel(RF),h,['RF series #' int2str(k)])
        i_min = 8*(k - 1);
        subaper=(1:32)+i_min;

        txdel_comp = NaN(1,128);
        txdel_comp(subaper) = txdel_suba; % Délais pour shoot seulement 32 des 128 elems

        RF{k}=simus(xs,zs,RC,txdel_comp,param,option); 

        txdel{k}=txdel_comp;
    end
    %close(h)
    
    % IQ demodulation
    IQ = cell(numel(RF),1);  % this cell will contain the I/Q series
    
    for k = 1:numel(IQ)
        IQ{k} = rf2iq(RF{k},param.fs,param.fc);
    end
    
    
    % Formation des voies 
    param.fnumber = [];
    [xi,zi] = impolgrid([256 256],l_pict,param);
    
    bIQ = zeros(256,256,numel(IQ));  % this array will contain the 13 images
    
    %h = waitbar(0,'Démodulation DAS');
    for k = 1:numel(IQ)
        %waitbar(k/numel(IQ),h,['DAS: I/Q series #' int2str(k)])
        i_min = 8*(k - 1);
        valid_indices = (1:32)+i_min;
        %waitbar(k/numel(IQ),h,['DAS: I/Q series #' int2str(k)])
        bIQ(valid_indices,:,k) = das(IQ{k},xi(valid_indices,:),zi(valid_indices,:),txdel{k},param);
    end
    %close(h)

    % Time Gain compensation
    %bIQ = tgc(bIQ);
    
    % Formation de l'image en BMODE par somme des images
    bIQ_sum = sum(bIQ, 3);
    I = bmode(bIQ_sum(:,:,1),40); % log-compressed image

