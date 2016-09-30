clear;

D=24; 
%uu=-5; potential=10; interaction=0; 

uu = 1;
potential = -1;
interaction = 10;


%gap =0.1 available Ds 4 20 26 32 38
%D=20;
%potential =  sqrt(1+(0.1)^2);
%interaction = 0;
%uu = -1/2;

%D=18;
%potential =  -1;
%interaction = 10;
%uu = 1;

fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);
load(fIn);
l = reshape(eye(D), 1, D^2); 
r = reshape(transpose(matr), D^2, 1);

fOut = sprintf('data_LLuu/D=%d/excitation_data/cpx_excitationsD=%duu=%dv=%dc=%d.mat', D , D, uu, potential, interaction);

p_min = -10; p_max = 10;

delta_p=0.2;

momenta = p_min:delta_p:p_max;

%H = cell(1, numel(plot_x));
%Heigvectors = cell(1, numel(plot_x));

%calculateExcitations(p, R, Q, l, r, matr, Tminus_in, D, mass, interaction, potential, uu, FOR_PLOT)

[excitation_energies, H, Heigvectors] = calculateExcitations(momenta, R, Q, l, r, matr, 0, D, mass, interaction, potential, uu, 1);


%save(fOut);
