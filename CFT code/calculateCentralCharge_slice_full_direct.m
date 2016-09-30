%%calulate central charge directly by calculating entropy as a function of
%%position of a strip

clear;

D=12;

potential = 10;
interaction = 0;
uu = -5;

symGauge = false;

scale = 0.001:0.1:10;

fIn = sprintf('data_LLuu/lr_eigenvectors/uu=%dv=%dc=%deigenvalues_all.mat', uu, potential, interaction);

load(fIn, 'dim', 'eigvals');

entropy = zeros(1, numel(scale));

centralCharge = 1; %expected central charge...

kappa = 6/( centralCharge*(sqrt(12/centralCharge)+1) );



%do interpolation starting for all ranges starting from eigval_start:dim(end) to eigval_start+scan_end:dim(end) 
scanDelta = 10;

%eigvals is sorted w.r.t. abs, uncomment this to sort using real part
for pp=1:numel(dim)
   eigvals{pp}   = -1*sort(-1*real(eigvals{D-3}));
end

mu = eigvals{D-3}(2);
corrLength = real(-1/mu);
position = scale*corrLength;


    fIn =sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);
    load(fIn, 'R', 'Q', 'matr' ); matl = eye(D);
    if symGauge == true
        symmetricGauge;
    end


for pp=1:numel(scale)
    
      fprintf('%d ', scale(pp));
        if (rem(pp,10)==0)
                fprintf('\n');
        end 
  
    
    
    TT = kron(R, conj(R)) + kron(Q, eye(D)) + kron(eye(D), conj(Q));
    
    expTTtmp = expm((position(pp))*TT);
    
    clearvars TT;
    
    expTTtmp = reshape(transpose(expTTtmp), D, D, D, D);
    
    expTTtmp2 = expTTtmp;
    for kk=1:D
        for mm = 1:D
            expTTtmp2(:, kk, mm, :) = expTTtmp(:, mm,  kk, :);
        end
    end
    expTTtmp = expTTtmp2;
    
    clearvars expTTtmp2;
    
    expTT = reshape((expTTtmp), D^2, D^2);
    
    clearvars expTTtmp;
    
    %%%%%%%%

    tmpmat = kron((transpose(matl^(0.5))), (transpose(matr))^(0.5))*(transpose(expTT))*kron((matl)^0.5, transpose(matr)^0.5);
    schmidt = svd(tmpmat);
    
    clearvars expTT tmpmat;

    entropy(pp) = -sum(schmidt.*log(schmidt));
    
end

fprintf('\n');


figure; plot((scale), entropy, 'x');

scale_log = log(scale);

entropy_pp = interp1(scale_log, entropy, 'spline','pp');

c_estimate_pp =fnder(entropy_pp,1);


c_estimate=3*ppval(c_estimate_pp, scale_log);

figure; plot(scale, c_estimate, 'x');


