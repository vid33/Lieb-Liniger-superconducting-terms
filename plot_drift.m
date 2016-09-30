clear;

power = 6;

E_eff = [ ] ;

c_eff = [ ] ;

D =  [  ];


for (m =1:power)
    
    fIn = sprintf('data/D=%d/cpxD=%dv=%dc=%d.mat', 2^m, 2^m, 1, 1);
    load(fIn, 'interaction_effective', 'energyDensity_effective');
    
    E_eff = [ E_eff, energyDensity_effective];
    c_eff = [ c_eff, interaction_effective];
    
    D = [ D 2^m ] ;
end

plot(D(1:end), E_eff(1:end), 'xk');

figure; plot(D(1:end), c_eff(1:end), 'xk');

