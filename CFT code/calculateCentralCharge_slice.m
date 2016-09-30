clear;

potential = 10;
interaction = 0;
uu = -5;

symGauge = true;

dimStart = 4;  dimEnd = 6;

dim = dimStart:1:dimEnd;

centralCharge = 1; %expected central charge...

kappa = 6/( centralCharge*(sqrt(12/centralCharge)+1) );

scale  =1;

%eigenvectorMax = 50;

%l2 = cell(eigenvectorMax, numel(dim));
%r2 = cell(eigenvectorMax, numel(dim));

%matl2 = cell(eigenvectorMax, numel(dim));
%matr2 = cell(eigenvectorMax, numel(dim));
%mu = cell(eigenvectorMax, numel(dim));

expTT = cell(1, numel(dim));
expTT2 = cell(1, numel(dim));
eigTT = cell(1, numel(dim));

entropy = zeros(1,  numel(dim));
entropy0 = zeros(1,  numel(dim));
%entropy2 = zeros(1,  numel(dim));
schmidt = cell(1,  numel(dim));
schmidt2 = cell(1,  numel(dim));

mu2 = zeros(1, numel(dim));

for pp=1:numel(dim)
    
      fprintf('%d ', dim(pp));
        if (rem(dim(pp),10)==0)
                fprintf('\n');
        end 
    
    fIn =sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', dim(pp), dim(pp), uu, potential, interaction);
    load(fIn);
    if symGauge == true
        symmetricGauge;
    end
    
    D = dim(pp);

    TT = kron(R, conj(R)) + kron(Q, eye(D)) + kron(eye(D), conj(Q));
    eigTT{pp} = sort(eig(TT));
    mu2(pp) = eigTT{pp}(3);
    
    expTTtmp = expm((-scale/mu2(pp))*TT);
    expTT2{pp} = expTTtmp;
    
    expTTtmp = reshape(transpose(expTTtmp), D, D, D, D);
    
    expTTtmp2 = expTTtmp;
    for kk=1:D
        for mm = 1:D
            expTTtmp2(:, kk, mm, :) = expTTtmp(:, mm,  kk, :);
        end
    end
    expTTtmp = expTTtmp2;
    
    expTTtmp = reshape((expTTtmp), D^2, D^2);
    
    expTT{pp} = expTTtmp;

end


fprintf('\n');

for pp=1:numel(dim)
    
    D = dim(pp);
    
    fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);
    load(fIn); matl = eye(D);
   
    if symGauge == true
        symmetricGauge;
    end
    
    
      fprintf('%d ', dim(pp));
        if (rem(dim(pp),10)==0)
                fprintf('\n');
        end 

   % densityMatrix = zeros(D^2); 
   % for kk=1:eigenvectorMax
   %     densityMatrix = densityMatrix + exp((mu{2,pp}))*kron(transpose(matl2{kk, pp}), matr2{kk, pp});
   % end

   %tmpmat = kron((transpose(matl^(0.5))), (transpose(matr))^(0.5))*(((transpose(expTT{pp}))^0.5));
   tmpmat = kron((transpose(matl^(0.5))), (transpose(matr))^(0.5))*(transpose(expTT{pp}))*kron((matl)^0.5, transpose(matr)^0.5);
   %tmpmat = tmpmat*tmpmat;
    schmidt{pp} = svd(tmpmat);
    
    %Z1 = trace(expTT2{pp});
    %schmidt2{pp} = svd(sqrt(expTT2{pp})/Z1);

    entropy(pp) = -sum(schmidt{pp}.*log(schmidt{pp}));
   % entropy2(pp) = -sum(schmidt2{pp}.*log(schmidt2{pp}));
    

end

%figure; plot(log(1./mu2), entropy0, 'x');
figure; plot(log(-1./mu2), entropy, 'xb-', log(-1./mu2), 2*entropy0, 'xr-');

%figure; plot(log(dim), entropy, 'xb-', log(dim), 2*entropy0, 'xr-');

