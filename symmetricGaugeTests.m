potential = 10;
uu=-5;
interaction = 0;

dim = (4:16);

TT_diag = zeros(1, numel(dim));

eig_TT_before = cell(1, numel(dim));
eig_TT_after = cell(1, numel(dim));

for kk=1:numel(dim)
    
    D = dim(kk);
    
    fprintf('%d ', kk);
    if (rem(kk,10)==0)
                fprintf('\n');
    end 
   
    fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', dim(kk), dim(kk), uu, potential, interaction); 
    load(fIn, 'energyDensity_new','matr', 'R', 'Q', 'Rkin','normRtangent');
    
    matl = eye(D);
    
    TT_before = kron(R, conj(R)) + kron(Q, eye(D)) + kron(eye(D), conj(Q));
    
    eig_TT_before{kk}  = eig(TT_before);
    
    symmetricGauge;
    TT = kron(R, conj(R)) + kron(Q, eye(D)) + kron(eye(D), conj(Q));
    
    eig_TT_after{kk}  = eig(TT);
    
    TT_diag_tmp = sort(diag(TT));

    TT_diag(kk) = TT_diag_tmp(1);
    
end