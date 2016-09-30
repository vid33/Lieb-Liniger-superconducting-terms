uu = 1;
potential = -1;
interaction = 10;

kappa = 6/(sqrt(12) + 1);

dim = 2:25;

%dim = 2;

massGap = zeros(1, numel(dim));

secondExcitation = zeros(1, numel(dim));

corrLength = zeros(1, numel(dim));

momentum = 0;

for kk=1:numel(dim)
    
    fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', dim(kk), dim(kk), uu, potential, interaction);
    load(fIn, 'D', 'R', 'Q', 'l', 'r', 'matr', 'mass');
    
    
%    [excitation_energies, H, Heigvectors] = calculateExcitations(p, R, Q, l, r, matr, Tminus_in, D, mass, interaction, potential, uu, FOR_PLOT)
   
    [excitation_energies, H, Heigvectors] = calculateExcitations(momentum, R, Q, l, r, matr, 0, D, mass, interaction, potential, uu, 0);
    
    excitation_energies = sort(excitation_energies);
    
    TT = kron(R, conj(R)) + kron(Q, eye(D)) + kron(eye(D), conj(Q));
    
    TT_eig = sort(eig(TT));
    
    corrLength(kk) = -1/TT_eig(2);
    
    massGap(kk) = min(excitation_energies);
    
    secondExcitation(kk) = excitation_energies(2);
    
    %[excitation_energies, H, Heigvectors] = calculateExcitations_phase(p, R, Q, l, r, matr, Tminus_in, D, mass, interaction, potential, uu, phase, FOR_PLOT)
    
    %[excitation_energies, H, Heigvectors] = calculateExcitations_phase(-pi/occulationDensity, R, Q, l, r, matr, 0, D, mass, interaction, potential, uu, -1, 0);
    
   % excitation_energies = sort(excitation_energies);
    
   % massGap_phase(kk) = min(excitation_energies);
    
end

figureName = sprintf('Mass Gap vs. D');
figure('Name', figureName);
        xlabel('D'); 
        ylabel('mass gap');
        hold all;
plot(dim, massGap, '-');
plot(dim, secondExcitation, '-');
hold off;

figure; plot(dim, corrLength, 'x');

%DATA FOR TOPOLOGICALLY NON-TRIVIAL:
%uu = 1; potential = -1; interaction = 10;
dim2 = [2, 3, 4, 5, 6, 7, 8, 10 ,12, 14, 16, 18];

massGap_phase = [0.47524, 0.47907, 0.486198, 0.486056, 0.486033, 0.485951, 0.485964, 0.486966, 0.485964, 0.485964, 0.485963, 0.485962];

hold all;
%plot(dim, 2*massGap_phase, 'x');
plot(dim2, 2*massGap_phase, '-');
hold off;
        xlabel('D'); 
        ylabel('mass gap');
        
       
