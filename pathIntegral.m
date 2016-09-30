clear;

D=8;

potential = 10;
interaction = 0;
uu = -5;

fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);

load(fIn);

J = 0.0000001;

TT = kron(R, conj(R)) + kron(eye(D), conj(Q)) + kron(Q, eye(D));

TT_J = (1 + J)*kron(R, conj(R)) + kron(eye(D), conj(Q)) + kron(Q, eye(D));

%TT_J = kron(R, conj(R)) + kron(eye(D), conj(Q)) + kron(Q, eye(D)) + J*(kron(R, eye(D)) + kron(eye(D), conj(R)));

masses = real(sort(eig(TT)));
masses_J = real(sort(eig(TT_J)));

%masses = masses.^2; masses_J = masses_J.^2;
masses = abs(masses); masses_J = abs(masses_J);

masses_delta = (masses - masses_J);

masses_der = masses_delta/J;
