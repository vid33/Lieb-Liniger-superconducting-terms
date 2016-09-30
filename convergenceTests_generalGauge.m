clear;

%m^2 = 1 lambda = sqrt(20) --> uu=-10/4 potential = 30/4


%potential =  (30/4);
%interaction = 0;
%uu = -(10/4);

mass = 1/2;
potential = 1;
interaction = 1;
uu = 1;

%potential = -1;
%interaction = 10;
%uu = 1;

%dim = 32:64;
%dim = (4:64);

dim = [6 8 12 16 20 24 32 36 ];

energy = zeros(1, numel(dim));
condensate = zeros(1, numel(dim));
corrLength = zeros(1, numel(dim));

entropy = zeros(1, numel(dim));

normTangentVec = zeros(1, numel(dim));

minSchmidt = zeros(1, numel(dim));

for kk=1:numel(dim)
    
    fprintf('%d ', kk);
    if (rem(kk,10)==0)
                fprintf('\n');
    end 
   
    fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d_generalGauge.mat', dim(kk), dim(kk), uu, potential, interaction); 
    load(fIn, 'energyDensity_new','matr', 'matl', 'R', 'Q', 'Rkin','normRtangent');
    
    
    TT = kron(R, conj(R)) + kron(Q, eye(dim(kk))) + kron(eye(dim(kk)), conj(Q));
    aa = sort(eig(TT));
    
    corrLength(kk) = aa(2);
    
    condensate(kk) = trace(matl*R*matr);
    %e_delta = 2*trace(Rkin*matr*Rkin');
    %energy(kk) = (1/sqrt(20))*( real(energyDensity_new + e_delta) ); THIS
    %IS ZERO!! WHY?
    
    e_delta = 0;
    energy(kk) = ( real(energyDensity_new + e_delta) );
    normTangentVec(kk) = real(normRtangent);
    
    schmidt_sq = diag(matl).^2;
    entropy(kk) = -sum(schmidt_sq.*log(schmidt_sq));
    
    minSchmidt(kk) = matl(end, end);
    
end

figure; plot(dim, energy, 'x');
        xlabel('$D$', 'Interpreter', 'LaTex'); 
        ylabel('energy');
        
        
figure; plot(dim, real(condensate),  'x');
        xlabel('$D$', 'Interpreter', 'LaTex'); 
        ylabel('$tr(R*r)$', 'Interpreter', 'LaTex');  
        


figure; plot(dim, real(corrLength), 'x');
        xlabel('$D$', 'Interpreter', 'LaTex'); 
        ylabel('Correlation Length');   
        
%figure; plot(dim, normTangentVec, 'x');
%        xlabel('$D$', 'Interpreter', 'LaTex'); 
%        ylabel('norm Rtangent');   

figure; plot(dim, real(entropy), 'x');
        xlabel('$D$', 'Interpreter', 'LaTex'); 
        ylabel('Entropy'); 
        
figure; plot(dim, real(minSchmidt), 'x');
        xlabel('$D$', 'Interpreter', 'LaTex'); 
        ylabel('min. Schmidt'); 
        

