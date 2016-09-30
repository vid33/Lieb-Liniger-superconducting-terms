clear;

D=8;

potential = 1;
interaction = 1;
%uu_list = [-0.5 -0.2 -0.1 -0.01 0 0.01 0.1 0.2 0.3 0.4 0.5];

%uu_list = [-0.1 0 0.1 0.2 0.3 0.4 0.5 1 1.5 2 5];

uu_list = [0.8 0.85 0.9 1 1.09 1.099 1.1 1.105 1.106 1.107 1.108 1.109 1.11 1.12 1.15 1.2] ;

corrLength = zeros(1, numel(uu_list));
condensateNorm = zeros(1, numel(uu_list));

excitationEnergy_min = zeros(1, numel(uu_list));

condensateReal = zeros(1, numel(uu_list));
condensateImag = zeros(1, numel(uu_list));

for kk=1:numel(uu_list)
   
    fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu_list(kk), potential, interaction);

    load(fIn, 'R', 'Q', 'matr');
    
    TT = kron(R, conj(R)) + kron(Q, eye(D)) + kron(eye(D), conj(Q));
    matl = eye(D); r = reshape(transpose(matr), D^2,1); l = reshape(transpose(matl), 1, D^2);
    
    eigvals = sort(eig(TT));
    
    corrLength(kk) = real(-1/eigvals(2));

    condensateNorm(kk) = norm(trace(R*matr));
    
    %[excitation_energies, H, Heigvectors] = calculateExcitations(p, R, Q, l, r, matr, Tminus_in, D, mass, interaction, potential, uu, FOR_PLOT)
    Tminus = -kron(Q, eye(D))- kron( eye(D), conj(Q)) - kron(R, conj(R)) - kron(r,l) - 1i*0*eye(D^2);
    [excitation_energies, ~, ~] = calculateExcitations(0, R, Q, l, r, matr, Tminus, D, 1/2, interaction, potential, uu_list(kk), 0);
    
    excitationEnergy_min(kk) = min(excitation_energies);
    
 %   condensateReal(kk) = real(trace(R*matr));
 %   condensateImag(kk) = imag(trace(R*matr));
    
end

figureName = sprintf('corrLength vs. uu');
figure('Name', figureName); plot(uu_list, corrLength, 'x-');

figureName = sprintf('condensateNorm vs. uu');
figure('Name', figureName); plot(uu_list, condensateNorm, 'x-');

figureName = sprintf('min excitation energy vs. uu');
figure('Name', figureName); plot(uu_list, excitationEnergy_min, 'x-');

%figureName = sprintf('condensateReal vs. uu');
%figure('Name', figureName); plot(uu_list, condensateReal, 'x-');

%figureName = sprintf('condensateImag vs. uu');
%figure('Name', figureName); plot(uu_list, condensateImag, 'x-');