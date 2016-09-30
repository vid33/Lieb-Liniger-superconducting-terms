clear;

uu = -5;
potential = 10;
interaction = 0;

D=4;

eigenvectorNo_max = D^2;


mu = zeros(D^2, 1);

l2 = cell(D^2, 1);
r2 = cell(D^2, 1);

matl2 = cell(D^2, 1);
matr2 = cell(D^2,1 );

    
    fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);
 
    load(fIn, 'R', 'Q', 'D', 'matr', 'matl', 'r', 'l', 'energyDensity_new'); matl = eye(D); mass = 1/2;
    symmetricGauge;
    r = reshape(transpose(matr), D^2, 1);
        
    TT = kron(Q, eye(D)) + kron(eye(D), conj(Q)) + kron(R, conj(R)); 
    
    [Vr, Dr] = eig(TT);
    [Vl, Dl] = eig(TT.'); Vl = conj(Vl);
    
    [Dr, Dr_index] = sort(diag(Dr));
    [Dl, Dl_index] = sort(diag(Dl));
    
    TTtest = zeros(D^2, D^2);
    
    for zz=1:D^2
     
        r2{zz} = Vr(:, Dr_index(zz));
        matr2{zz} = transpose(reshape(r2{zz}, D, D));
        
        l2{zz} = (Vl(:, Dl_index(zz)) )';
        matl2{zz} = transpose(reshape(l2{zz}, D, D));
        
        tmpNorm = trace(transpose(matl2{zz})*matr2{zz});
        
        matr2{zz} = matr2{zz}/tmpNorm;
        r2{zz} = r2{zz}/tmpNorm;
        
        mu(zz) = Dl(zz);
        
        TTtest = TTtest + mu(zz)*kron(r2{zz}, l2{zz});
        
    end  
    
 %   [ l2{k}, r2{k}, matl2{k}, matr2{k}, mu(k), r2_eigval_test] = lrEigenvector_subleading( R, Q, D, eigenvector);
    
    
scale = 2;
 x_max = real(-1/mu(2))*scale;
 x_min = 0.1;
plot_x = linspace(x_min, x_max, 50);
correlator = zeros(1, numel(plot_x));

for kk = 1:numel(plot_x)
correlator(kk) = l*(kron(R, eye(D)) - kron(eye(D), conj(R)) )*expm(plot_x(kk)*TT)*( kron(R, eye(D)) -kron(eye(D), conj(R)))*r;
end

testOverlap1 = zeros(1, numel(D^2));
testOverlap2 = zeros(1, numel(D^2));
OPcoeffs = zeros(numel(D^2), numel(D^2));

OPtest = zeros(D^2, D^2);
TTtest2 = zeros(D^2, D^2);

OP = kron(R, eye(D)) - kron(eye(D), conj(R));
Rkin = Q*R - R*Q;

%OP = (1/(2*mass))*kron(Rkin, conj(Rkin)) ...
%    + uu*(kron(R*R, eye(D)) + kron(eye(D), conj(R)*conj(R))) ...
%    + potential*kron(R, conj(R));

for kk=1:D^2
testOverlap1(kk) = l2{kk}*OP*r;
testOverlap2(kk) = l*OP*r2{kk};

end

for kk=1:D^2
    for mm=1:D^2
            l2{kk}*r2{mm};
          OPcoeffs(kk,mm) = l2{kk}*OP*r2{mm};
           %l2{kk}*TT*r2{mm}
           %TTtest2 = TTtest2+l2{kk}*TT*r2{mm}*kron(r2{mm}, l2{kk});
    end
end

for kk=1:D^2
    for mm=1:D^2

            OPtest = OPtest + OPcoeffs(kk,mm)*kron(r2{kk}, l2{mm});
    end
end


%figure; plot(1:1:D^2, testOverlap1, 'x', 1:1:D^2, testOverlap2, 'x');

testOverlap_sq = testOverlap1.*testOverlap2;

figure; plot(1:1:D^2, (testOverlap_sq),'x');

%figure; plot(1:1:D^2, testOverlapDensity, 'x');




