   clear;

   D = 8;
   
    uu = -5;
    potential = 10;
    interaction = 0;

    fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);
 
    load(fIn);
    
    schmidt = svd(matr);
    entropy = sum(-schmidt.*log(schmidt));
    
    
    TT = kron(R, conj(R)) + kron(Q, eye(D)) + kron(eye(D), conj(Q));
    
    %symmetricGauge;
    
    schmidtTT = svd(expm(1000*TT));
    
    entropyTT =  sum(-(schmidtTT(1)^2)*log((schmidtTT(1)^2)));
    
    