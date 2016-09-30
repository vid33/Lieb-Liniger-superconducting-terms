%D=24;

%uu = 1;
%potential = -1;
%interaction = 10;

D=16; uu=-5; potential = 10; interaction=0;


    fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);
    load(fIn);


TT = kron(R, conj(R))  +  kron(Q, eye(D)) + kron(eye(D), conj(Q));
 
 eigvals = sort(eig(TT));
 
 figure; plot(imag(eigvals), -real(eigvals), 'x');