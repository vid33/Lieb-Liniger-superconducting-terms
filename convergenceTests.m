clear;

%m^2 = 1 lambda = sqrt(20) --> uu=-10/4 potential = 30/4


%potential =  (30/4);
%interaction = 0;
%uu = -(10/4);

mass = 1/2;
potential = 1;
interaction = 1;
uu = 3;

%potential = -1;
%interaction = 10;
%uu = 1;

%dim = 32:64;
%dim = (4:64);

dim = [2 4 6 8 10 12 14 16 18 20 24];

energy = zeros(1, numel(dim));
condensate = zeros(1, numel(dim));
corrLength = zeros(1, numel(dim));

RtangentNorm = zeros(1, numel(dim));

for kk=1:numel(dim)
    
    fprintf('%d ', kk);
    if (rem(kk,10)==0)
                fprintf('\n');
    end 
   
    fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', dim(kk), dim(kk), uu, potential, interaction); 
    load(fIn, 'energyDensity_new','matr', 'R', 'Q', 'Rkin','normRtangent');
    RtangentNorm(kk) = sqrt(real(normRtangent));

    
   TT = kron(R, conj(R)) + kron(Q, eye(dim(kk))) + kron(eye(dim(kk)), conj(Q));
   aa = sort(eig(TT));
    
   corrLength(kk) = -1/aa(2);
    
    condensate(kk) = trace(R*matr);
    %e_delta = 2*trace(Rkin*matr*Rkin');
    %energy(kk) = (1/sqrt(20))*( real(energyDensity_new + e_delta) ); THIS
    %IS ZERO!! WHY?
    
    e_delta = 0;
    energy(kk) = ( real(energyDensity_new + e_delta) );
    
end

figure; plot(dim, energy, 'x');
xlabel('$D$', 'Interpreter', 'LaTex'); 
ylabel('energy');


figure; plot(dim, real(condensate),  'x');
xlabel('$D$', 'Interpreter', 'LaTex'); 
ylabel('$tr(R*r)$', 'Interpreter', 'LaTex');  

figure; plot(dim, corrLength, 'x');
xlabel('$D$', 'Interpreter', 'LaTex'); 
ylabel('Correlation Length', 'Interpreter', 'LaTex'); 

figure; plot(dim, log10(RtangentNorm), 'x');
xlabel('$D$', 'Interpreter', 'LaTex'); 
ylabel('log10(norm of tangent vector)', 'Interpreter', 'LaTex');
        
        
       break; 
        
        
fIn= sprintf('data_LLuu/lr_eigenvectors/uu=%dv=%dc=%deigenvector=%d.mat',uu,  potential, interaction , 2);
load(fIn, 'mu');

figure; plot(dim, real(-1./mu), 'x');
        xlabel('$D$', 'Interpreter', 'LaTex'); 
        ylabel('Correlation Length');   
        
        figure; plot(dim, normTangentVec, 'x');
        xlabel('$D$', 'Interpreter', 'LaTex'); 
        ylabel('norm Rtangent');   
        

