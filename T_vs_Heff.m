clear;

kappa = 6/(sqrt(12) + 1);
%kappa = 1.2945;

Lambda = sqrt(20);

D=16; uu=-5; potential = 10; interaction=0;

fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D , D, uu, potential, interaction);
load(fIn); occupationDensity = real(occupationDensity);

TT = kron(Q, eye(D)) + kron(eye(D), conj(Q)) + kron(R, conj(R));

TT_eigvals = real(-1*sort(eig(TT)));

k_from_TT = sqrt(sqrt( TT_eigvals.^2 + Lambda^4/4) - Lambda^2/2);

E_zero_at_TT = k_from_TT/Lambda;

excitation_TT = TT_eigvals - E_zero_at_TT;


k_exact  = -25:0.1:25;
omega_cutoff = sqrt( (k_exact.^2 + Lambda^2/2).^2 - (1/4)*Lambda^4);

omega_exact = abs(Lambda*k_exact);

figure; plot(k_exact, omega_cutoff);
hold all;
plot(k_exact, omega_exact);

hold off;

%l = reshape(transpose(eye(D)), 1, D^2);
%r =reshape(transpose(matr), D^2, 1);

p = 4;
TT_p = TT + p*kron(r, l);

[excitation_energies, H, Heigvectors] = calculateExcitations(0, R, Q, l, r, matr, 0, D, mass, interaction, potential, uu, 0);

%calculateExcitations_phase(p, R, Q, l, r, matr, Tminus_in, D, mass, interaction, potential, uu, phase, FOR_PLOT)

%[excitation_energies, H, Heigvectors] = calculateExcitations_phase(0, R, Q, l, r, matr, 0, D, mass, interaction, potential, uu, -1, 0);



dim_squared = 0:D^2-1;
 figure; plot(dim_squared(1:end-1), excitation_energies(2:end)-excitation_energies(1), 'x');
 
 hold all;
 
 plot(dim_squared(1:end-1), TT_eigvals(2:end), 'o');
 
 hold off;

 
 figure; plot(log(TT_eigvals(2:end)), log(excitation_energies(2:end)-excitation_energies(1)), 'x');
 
  %figure; plot(log(TT_eigvals(2:300)), log(excitation_energies(2:300)-excitation_energies(1)), 'x');
 
 
