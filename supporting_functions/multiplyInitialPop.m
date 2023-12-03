function mat_sats = multiplyInitialPop(cfgMC,mat_sats,g1,g2,g3)
    inmatsat = mat_sats;
    mult = cfgMC.initpopMultiplier;
    % randomly choose which satellites to add to initial population (per object type)
    % note that physical parameters will be filled in by resampling (fillMassRadiusResample)
    if mult < 1  % simply select a subset of original matsats
        randind1 = randi(numel(g1.allclass), [floor(numel(g1.allclass) * mult),1]); % repeats allowed
        randind2 = randi(numel(g2.allclass), [floor(numel(g2.allclass) * mult),1]);
        randind3 = randi(numel(g3.allclass), [floor(numel(g3.allclass) * mult),1]);
        mat_sats = inmatsat([g1.allclass(randind1);g2.allclass(randind2);g3.allclass(randind3)],:); 
    elseif mult > 1 % if mult>1, ensure all of original matsats are included
        randind1 = randi(numel(g1.allclass), [floor(numel(g1.allclass) * (mult-1)),1]);
        randind2 = randi(numel(g2.allclass), [floor(numel(g2.allclass) * (mult-1)),1]);
        randind3 = randi(numel(g3.allclass), [floor(numel(g3.allclass) * (mult-1)),1]);
        extrasats = inmatsat([g1.allclass(randind1);g2.allclass(randind2);g3.allclass(randind3)],:); 
        % zero out mass and radius and fill in (resample)
        idx_mass = 8; idx_radius = 9;
        extrasats(:,[idx_mass, idx_radius]) = 0; 
        if cfgMC.fillMassRadius > 0
            addmatsats = fillMassRadiusResample(extrasats,g1,g2,g3);
        else
            addmatsats = extrasats;
        end
        % scramble argperigee and mean motion
        idx_argpo = 5; idx_mo = 6;  % range: [0,2pi]
        addmatsats(:,[idx_argpo, idx_mo]) = 2 * pi * rand(size(addmatsats,1),2);
        mat_sats = [mat_sats; addmatsats]; 
    end
end
