clear;

kappa = 6/(sqrt(12) + 1);
%kappa = 1.2945;

Lambda = sqrt(20);

dim = 10:32;

excitation_energies_all = cell(1, numel(dim));

excitation1 = zeros(1, numel(dim));
excitation2 = zeros(1, numel(dim));

for zz=1:numel(dim)
    
D=dim(zz); uu=-5; potential = 10; interaction=0;

fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D , D, uu, potential, interaction);
load(fIn); occupationDensity = real(occupationDensity);

TT = kron(R, conj(R)) + kron(eye(D), conj(Q)) + kron(Q, eye(D));

TT_eigvals = -1*sort(eig(TT));



[excitation_energies, H, Heigvectors] = calculateExcitations(0, R, Q, l, r, matr, 0, D, mass, interaction, potential, uu, 0);

excitation1(zz) = excitation_energies(2);
excitation2(zz) = excitation_energies(3);

excitationT1(zz) = TT_eigvals(2);
excitationT2(zz) = TT_eigvals(3);

excitation_energies_all{zz} = excitation_energies;
%calculateExcitations_phase(p, R, Q, l, r, matr, Tminus_in, D, mass, interaction, potential, uu, phase, FOR_PLOT)

%[excitation_energies, H, Heigvectors] = calculateExcitations_phase(0, R, Q, l, r, matr, 0, D, mass, interaction, potential, uu, -1, 0);

end


figure; plot(dim, excitation2./excitation1, 'x');

figure; plot(dim, excitationT2./excitationT1, 'x');

