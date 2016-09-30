clear;

D=2;

potential = 10;
interaction = 0;
uu = -5;

fIn =sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);
load(fIn);

symmetricGauge;
    
beta = 0.1:0.01:3;
alpha = 0.1;
    
ZZho = ones(1, numel(beta));
ZZsimple = zeros(1, numel(beta));

ZZtest = zeros(1, numel(beta));

    TT = kron(R, conj(R)) + kron(eye(D), conj(Q)) + kron(Q, eye(D));
    
    TTalpha = TT+ alpha*kron(R, eye(D)) + alpha*kron(eye(D), conj(R));
    
    TTeigs = -1*real(sort(eig(TT)));
    TTeigs_alpha = -1*real(sort(eig(TTalpha)));

for kk=1:numel(beta)
    
    ZZtest(kk) = trace(expm(beta(kk)*TT))-1;

    
    for mm=1:(numel(TTeigs)-1)
       
        ZZsimple(kk) = ZZsimple(kk) + exp(-beta(kk)*TTeigs(mm+1));
        
        ZZho(kk) = ZZho(kk)*exp(-beta(kk)*TTeigs(mm+1)/2)/(1- exp(-beta(kk)*TTeigs(mm+1))) ;
        
    end
    
    
end 

figure; plot(beta, ZZsimple, beta, ZZho);
    