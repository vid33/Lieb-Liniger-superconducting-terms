clear;

%D=16; uu=-5; potential=10; interaction=0; 

%gap =0.1 available Ds 4 20 26 32 38
%D=4;
%potential =  sqrt(1+(0.1)^2);
%interaction = 0;
%uu = -1/2;

D=16;
potential =  -1;
interaction = 10;
uu = 1;

phase = -1;

fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);
load(fIn);

fOut = sprintf('data_LLuu/D=%d/excitation_data_phase/cpx_excitationsD=%duu=%dv=%dc=%dphase=%d.mat', D , D, uu, potential, interaction, phase);

p_min = -3; p_max = 3;

delta_p=0.1;

momenta = p_min:delta_p:p_max;

%H = cell(1, numel(plot_x));
%Heigvectors = cell(1, numel(plot_x));

%calculateExcitations(p, R, Q, l, r, matr, Tminus_in, D, mass, interaction, potential, uu, FOR_PLOT)

%[excitation_energies, H, Heigvectors] = calculateExcitations(momenta, R, Q, l, r, matr, 0, D, mass, interaction, potential, uu, 1);

[excitation_energies, H, Heigvectors] = calculateExcitations_phase(momenta, R, Q, l, r, matr, 0, D, mass, interaction, potential, uu, phase, 1);


%[excitation_energies, H, Heigvectors] = calculateExcitations_phase(p, R, Q, l, r, matr, Tminus_in, D, mass, interaction, potential, uu, phase, FOR_PLOT)


%save(fOut);