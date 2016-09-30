clear;

potential = 10;
interaction = 0;
uu = -5;

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
    
%    calculateFR( R, Rkin, Q, matl, matr, potential, uu, interaction, F_zero, D )

    [ FR, flag ] = calculateFR( R, Rkin, Q, matl, matr, potential, uu, interaction, Fzero, D );
    matFR = transpose(reshape(FR, D, D));

    [VmatFR, DmatFR] = eig(matFR, matr);
    
    [DmatFR, DmatFR_index] = sort(diag(DmatFR));
    
    %for zz=1:D
     
    VLsqrt = VmatFR(:, DmatFR_index(1));

    VLsqrt  = transpose(VLsqrt);

    matVL  = kron(VLsqrt', VLsqrt);

    VL = reshape(transpose(matVL), 1, D^2);
    
    norm = trace(transpose(matVL)*matr);
    matVL = matVL/norm;
    VL = VL/norm;
    VLsqrt = VLsqrt/sqrt(norm);
    
    
   fOut =  sprintf('data_LLuu/D=%d/boundary/boundaryL_D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);
    
   save(fOut, 'VL', 'matVL', 'VLsqrt', 'DmatFR');
   
   clearvars VL matVL VLsqrt R Q Rkin matr F_zero;
   
end

fprintf('\n');
