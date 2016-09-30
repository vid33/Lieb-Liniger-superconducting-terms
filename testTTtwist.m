clear;

D=4;

potential = 10;
interaction = 0;
uu = -5;

fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);

load(fIn);

%symmetricGauge;

TT = kron(R, conj(R)) + kron(Q, eye(D)) + kron(eye(D) , conj(Q));


Rcovec = reshape(transpose(R),D^2, 1);
RbarVec = reshape(transpose(conj(R)), 1, D^2);
Qcovec = reshape(transpose(Q), D^2, 1);
QbarVec = reshape(transpose(conj(Q)), 1, D^2);


TTtwist = kron(Rcovec, RbarVec) + kron(Qcovec, reshape(eye(D), 1, D^2)) + kron(reshape(eye(D), D^2, 1), QbarVec);

