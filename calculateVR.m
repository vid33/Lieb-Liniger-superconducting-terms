clear;

potential = 10;
interaction = 0;
uu = -5;

mass = 0.5;

precision = 1e-14;

dim = 2:1:64;

for zz = 1: numel(dim)
    
    fprintf('%d ', dim(zz));
    if (rem(zz,10)==0)
            fprintf('\n');
    end 
    
    D = dim(zz);
    
     fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);

    load(fIn, 'R', 'Q', 'matr'); matl = eye(D);
    Rkin = Q*R - R*Q;
    Fzero = cpxrand(D^2, 1);
    
    %function [ F, matF ] = calculateF( R, Rkin, Q, matr, potential, mass, uu, interaction, Fzero, D, precision )


    [ FL, matFL] = calculateF( R, Rkin, Q, matr, potential, mass, uu, interaction, Fzero, D, precision);
    %matFL = transpose(reshape(FL, D, D));

    [VmatFL, DmatFL] = eig(matFL);
    
    [DmatFL, DmatFL_index] = sort(diag(DmatFL));
    
    %for zz=1:D
     
    VRsqrt = VmatFL(:, DmatFL_index(1));

    %VRsqrt  = transpose(VRsqrt);

    matVR  = kron(VRsqrt', VRsqrt);
    matVR = conj(matVR);

    VR = reshape(transpose(matVR), D^2, 1);
    
    %norm = trace(transpose(matl)*matVR);
    %matVL = matVR/norm;
    %VR = VR/norm;
    %VRsqrt = VRsqrt/sqrt(norm);
    
    
   fOut =  sprintf('data_LLuu/D=%d/boundary/boundaryR_D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);
    
   save(fOut, 'VR', 'matVR', 'VRsqrt', 'DmatFL');
   
   clearvars VR matVR VRsqrt R Q Rkin matr F_zero;
   
end

fprintf('\n');
