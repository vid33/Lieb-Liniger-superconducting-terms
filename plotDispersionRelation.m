clear;

D=24;
%potential =  10;
%interaction = 0;
%uu = -5;

uu = 1;
potential = -1;
interaction = 10;


fIn = sprintf('data_LLuu/D=%d/excitation_data/cpx_excitationsD=%duu=%dv=%dc=%d.mat', D , D, uu, potential, interaction);
load(fIn);

l = reshape(eye(D), 1, D^2); 
r = reshape(transpose(matr), D^2, 1);

excitation_end = 3; 
momentum_end = numel(momenta);

plot_y = zeros(excitation_end, momentum_end);

for k = 1: momentum_end

    tmp = excitation_energies(1:D^2, k);
    
    tmp = sort(tmp);
    
    tmp = tmp(1:excitation_end);
    
    plot_y(:, k)  = tmp;
    
end

figureName = sprintf('D=%d, v=%g, c=%g', D, -potential, interaction);


figure('Name', figureName);  plot(real((1/(pi*occupationDensity))*momenta(1:momentum_end)), real( (1/((pi*occupationDensity)^2))*plot_y ),'-');

    xlabel('$p$', 'Interpreter', 'LaTex'); 
    ylabel('$\Delta E$', 'Interpreter', 'LaTex');

%plot(real(momenta(1:momentum_end)), real(plot_y ),'x');

%plot(plot_x(:,1 :25), plot_y(:, 1: 25), 'xk');

