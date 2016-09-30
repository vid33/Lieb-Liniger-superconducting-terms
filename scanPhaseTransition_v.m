clear;

D=8;

uu = 1;
interaction = 1;

v_list = [0.8 0.85 0.8684 0.87 0.871 0.8715 0.8717 0.8718 0.8719 0.872 0.875 0.88 0.9 1 1.1 ] ;

v_list_right = fliplr(0.87140:0.00001:0.87175);
v_list_left = [0.87151 0.87152 0.871525 0.87153  0.87154 ];

corrLength = zeros(1, numel(v_list));
corrLength_right = zeros(1, numel(v_list_right));
corrLength_left = zeros(1, numel(v_list_left));

condensateNorm = zeros(1, numel(v_list));
condensateNorm_right = zeros(1, numel(v_list_right));
condensateNorm_left = zeros(1, numel(v_list_left));

entropy = zeros(1, numel(v_list));
entropy_right = zeros(1, numel(v_list_right));
entropy_left = zeros(1, numel(v_list_left));

energy = zeros(1, numel(v_list));
energy_right = zeros(1, numel(v_list_right));
energy_left = zeros(1, numel(v_list_left));

excitationEnergy_min = zeros(1, numel(v_list));
excitationEnergy_min_left = zeros(1, numel(v_list_left));
excitationEnergy_min_right = zeros(1, numel(v_list_right));



for kk=1:numel(v_list)
   
    fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, v_list(kk), interaction);

    load(fIn, 'R', 'Q', 'matr', 'energyDensity_new');
    matl = eye(D); r = reshape(transpose(matr), D^2,1); l = reshape(transpose(matl), 1, D^2);
    
    TT = kron(R, conj(R)) + kron(Q, eye(D)) + kron(eye(D), conj(Q));
    
    eigvals = sort(eig(TT));
    
    corrLength(kk) = real(-1/eigvals(2));

    condensateNorm(kk) = norm(trace(R*matr));
    
     %[excitation_energies, H, Heigvectors] = calculateExcitations(p, R, Q, l, r, matr, Tminus_in, D, mass, interaction, potential, uu, FOR_PLOT)
    Tminus = -kron(Q, eye(D))- kron( eye(D), conj(Q)) - kron(R, conj(R)) - kron(r,l) - 1i*0*eye(D^2);
    [excitation_energies, ~, ~] = calculateExcitations(0, R, Q, l, r, matr, Tminus, D, 1/2, interaction, v_list(kk), uu, 0);
    
    excitationEnergy_min(kk) = min(excitation_energies);
    
    schmidt = svd(matr);
    entropy(kk) = sum(-schmidt.*log(schmidt));
    
    energy(kk) = real(energyDensity_new);
    
end

for kk=1:numel(v_list_right)
   
    fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d_right.mat', D, D, uu, v_list_right(kk), interaction);

    load(fIn, 'R', 'Q', 'matr', 'energyDensity_new');
    matl = eye(D); r = reshape(transpose(matr), D^2,1); l = reshape(transpose(matl), 1, D^2);
    
    TT = kron(R, conj(R)) + kron(Q, eye(D)) + kron(eye(D), conj(Q));
    
    eigvals = sort(eig(TT));
    
    corrLength_right(kk) = real(-1/eigvals(2));

    condensateNorm_right(kk) = norm(trace(R*matr));
    
    schmidt_right = svd(matr);
    entropy_right(kk) = sum(-schmidt_right.*log(schmidt_right));
    
    energy_right(kk) = real(energyDensity_new);
    
         %[excitation_energies, H, Heigvectors] = calculateExcitations(p, R, Q, l, r, matr, Tminus_in, D, mass, interaction, potential, uu, FOR_PLOT)
    Tminus = -kron(Q, eye(D))- kron( eye(D), conj(Q)) - kron(R, conj(R)) - kron(r,l) - 1i*0*eye(D^2);
    [excitation_energies, ~, ~] = calculateExcitations(0, R, Q, l, r, matr, Tminus, D, 1/2, interaction, v_list_right(kk), uu, 0);
    
    excitationEnergy_min_right(kk) = min(excitation_energies);
    
end

for kk=1:numel(v_list_left)
   
    fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d_left.mat', D, D, uu, v_list_left(kk), interaction);

    load(fIn, 'R', 'Q', 'matr', 'energyDensity_new');
    matl = eye(D); r = reshape(transpose(matr), D^2,1); l = reshape(transpose(matl), 1, D^2);
    
    TT = kron(R, conj(R)) + kron(Q, eye(D)) + kron(eye(D), conj(Q));
    
    eigvals = sort(eig(TT));
    
    corrLength_left(kk) = real(-1/eigvals(2));

    condensateNorm_left(kk) = norm(trace(R*matr));
    
    schmidt_left = svd(matr);
    entropy_left(kk) = sum(-schmidt_left.*log(schmidt_left));
    
    energy_left(kk) = real(energyDensity_new);
    
    %[excitation_energies, H, Heigvectors] = calculateExcitations(p, R, Q, l, r, matr, Tminus_in, D, mass, interaction, potential, uu, FOR_PLOT)
    Tminus = -kron(Q, eye(D))- kron( eye(D), conj(Q)) - kron(R, conj(R)) - kron(r,l) - 1i*0*eye(D^2);
    [excitation_energies, ~, ~] = calculateExcitations(0, R, Q, l, r, matr, Tminus, D, 1/2, interaction, v_list_left(kk), uu, 0);
    
    excitationEnergy_min_left(kk) = min(excitation_energies);

    
end

figureName = sprintf('corrLength vs. v');
figure('Name', figureName); plot(v_list, corrLength, 'x-');
hold on; plot(v_list_right, corrLength_right, 'xr-'); 
plot(v_list_left, corrLength_left, 'xg-');  hold off;

figureName = sprintf('condensateNorm vs. v');
figure('Name', figureName); plot(v_list, condensateNorm, 'x-');
hold on; plot(v_list_right, condensateNorm_right, 'xr-'); 
plot(v_list_left, condensateNorm_left, 'xg-');  hold off;

figureName = sprintf('entropy vs. v');
figure('Name', figureName); plot(v_list, entropy, 'x-');
hold on; plot(v_list_right, entropy_right, 'xr-'); 
plot(v_list_left, entropy_left, 'xg-'); hold off;

figureName = sprintf('energy vs. v');
figure('Name', figureName); plot(v_list, energy, 'x-');
hold on; plot(v_list_right, energy_right, 'xr-'); 
plot(v_list_left, energy_left, 'xg-'); hold off;

figureName = sprintf('min excitation energy vs. v');
figure('Name', figureName); plot(v_list, excitationEnergy_min, 'x-');
hold on; plot(v_list_right, excitationEnergy_min_right, 'xr-'); 
plot(v_list_left, excitationEnergy_min_left, 'xg-'); hold off;


