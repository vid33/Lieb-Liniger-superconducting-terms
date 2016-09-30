clear;

potential = 10;
interaction = 0;
uu = -5;
mass = 0.5;

dimStart = 5;  dimEnd = 64;

dim = dimStart:1:dimEnd;


energyL = zeros(1,  numel(dim));
energyR = zeros(1,  numel(dim));

ratio23 = zeros(1, numel(dim));

for pp=1:numel(dim)
    
      fprintf('%d ', dim(pp));
        if (rem(dim(pp),10)==0)
                fprintf('\n');
        end 
    
        D = dim(pp);
    
    fIn = sprintf('data_LLuu/D=%d/boundary/boundaryL_D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);
    load(fIn, 'DmatFR');
    fIn = sprintf('data_LLuu/D=%d/boundary/boundaryR_D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);
    load(fIn, 'DmatFL');
    
    energyL(pp) = DmatFR(1);
    energyR(pp) = DmatFL(1);
    
    ratio23(pp) = (DmatFR(4)-DmatFR(1))/(DmatFR(2) - DmatFR(1));
   
end

fprintf('\n');

figureName = sprintf('Energy vs. D');
figure('Name', figureName);
hold on;
plotL= plot(dim, energyL, 'xb');
plotR = plot(dim, energyR, 'xr');

legend([plotL,plotR],'energyL','energyR');

hold off;

figure; plot(dim, ratio23, 'xb');
