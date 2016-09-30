%expectation values in terms of derivatives wrt J

clear;

D=8;

potential = 10;
interaction = 0;
uu = -5;

fIn =sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);
load(fIn);

symmetricGauge;
    
beta = 0:0.001:3;
    
TT_beta_eig_max = zeros(1, numel(beta));

    TT = kron(R, conj(R)) + kron(eye(D), conj(Q)) + kron(Q, eye(D));
        

for kk=1:numel(beta)
    
    TT_beta = TT+ beta(kk)*( kron(R, eye(D)) + kron(eye(D), conj(R)) );
    
    eigvals = real(eig(TT_beta));
    
    eigvals_pos_index = find(eigvals > 0 );
    
    
    if numel(eigvals_pos_index) > 0
    
        eigvals_pos = zeros(1, numel(eigvals_pos_index));
    
        for zz = 1:numel(eigvals_pos_index)
            eigvals_pos(zz) = eigvals(eigvals_pos_index(zz));
        end
    
        eigvals_pos = sort(eigvals_pos);
    
        TT_beta_eig_max(kk) = eigvals_pos(end);
    else
        fprintf('Was HERE!\n');
        eigvals = sort(eigvals);
        TT_beta_eig_max(kk) = eigvals(end);
        
        %break;
    end
    
    
end 

figure; plot(beta, TT_beta_eig_max);
    