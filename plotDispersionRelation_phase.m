%clear;

%D=16; uu=-5; potential=10; interaction=0;

%gap =0.1 available Ds 4 20 26 32 38
%D=4;
%potential =  sqrt(1+(0.1)^2);
%interaction = 0;
%uu = -1/2;

D=16;
potential =  10;
interaction = 0;
uu = -5;

phase=-1;

fIn = sprintf('data_LLuu/D=%d/excitation_data_phase/cpx_excitationsD=%duu=%dv=%dc=%dphase=%d.mat', D , D, uu, potential, interaction, phase);


load(fIn);

excitation_end = 4; 
momentum_end = numel(momenta);

plot_y = zeros(excitation_end, momentum_end);

for k = 1: momentum_end

    tmp = excitation_energies(1:D^2, k);
    
    tmp = sort(tmp);
    
    tmp = tmp(1:excitation_end);
    
    plot_y(:, k)  = tmp;
    
end

figureName = sprintf('D=%d, uu-%g, v=%g, c=%g', D, uu, potential, interaction);

figure('Name', figureName);  
%plot(real((1/(occupationDensity))*momenta(1:momentum_end)) - pi, real( (1/((occupationDensity)^2))*plot_y ),'-');

plot(real(momenta(1:momentum_end)), real(plot_y ),'x');

%plot(plot_x(:,1 :25), plot_y(:, 1: 25), 'xk');

