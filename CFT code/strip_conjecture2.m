clear;

uu = -5;
potential = 10;
interaction = 0;


dim = 10;

centralCharge = 1; %expected central charge...

kappa = 6/( centralCharge*(sqrt(12/centralCharge)+1) );

scale  =0.5;

eigvals_tilde = cell(1, numel(dim));
ratio_orig = zeros(1, numel(dim));
ratio_tilde =zeros(1, numel(dim));
ratio_tilde_corrected =zeros(1, numel(dim));

mu2 = zeros(1, numel(dim));


fIn = sprintf('data_LLuu/lr_eigenvectors/uu=%dv=%dc=%deigenvalues_all.mat',uu, potential, interaction);

load(fIn, 'eigvals');

correlationEigenvalue = 2;
for pp=1:numel(dim)
  
    mu2(pp) = eigvals{dim(pp) -1}(correlationEigenvalue);
end
    
    
for pp=1:numel(dim)
    
      fprintf('%d ', dim(pp));
        if (rem(dim(pp),10)==0)
                fprintf('\n');
        end 
    
    fIn =sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', dim(pp), dim(pp), uu, potential, interaction);
    
    load(fIn);
    
    D = dim(pp);

   % TT = kron(R, conj(R)) + kron(Q, eye(D)) + kron(eye(D), conj(Q));
   % eigTT = sort(eig(TT));
   % mu2(pp) = eigTT(2);
  
    pos_in = real(-1/mu2(pp))*scale;
   
    TT = kron(R, conj(R)) + kron(Q, eye(D)) + kron(eye(D), conj(Q));
    
    expTTtmp = expm((-scale/mu2(pp))*TT);
    
    expTTtmp = reshape(transpose(expTTtmp), D, D, D, D);
    
    expTTtmp2 = expTTtmp;
    for kk=1:D
        for mm = 1:D
            expTTtmp2(:, kk, mm, :) = expTTtmp(:, mm,  kk, :);
        end
    end
    expTTtmp = expTTtmp2;
    
    expTT = reshape((expTTtmp), D^2, D^2);
    
    TT_tilde = logm(expTT);
    
    eigvals_tilde{pp} = sort(eig(TT_tilde));
    
    ratio_orig(pp) = real( eigvals{dim(pp)-1}(3)/eigvals{dim(pp)-1}(2) );
    
    ratio_tilde(pp) = real( eigvals_tilde{pp}(3)/eigvals_tilde{pp}(2) );
    
    ratio_tilde_corrected(pp) = real( (eigvals_tilde{pp}(3) - eigvals_tilde{pp}(1))/(eigvals_tilde{pp}(2) - eigvals_tilde{pp}(1)) );
    
end

figure; plot(dim, ratio_orig, '-');

figure; plot(dim, ratio_tilde, '-');
hold all;
plot(dim, ratio_tilde_corrected, '-');
hold off;
