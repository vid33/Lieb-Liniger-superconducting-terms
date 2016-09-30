clear;

uu = -5;
potential = 10;
interaction = 0;

dx = 0.0001;

D = 4;
fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);    
load(fIn); matl = eye(D);
%symmetricGauge; R = 0.5*(R+ R');
Rkin = Q*R -R*Q;

OP = R;

    
TT = kron(R, conj(R)) + kron(Q, eye(D)) + kron(eye(D), conj(Q));

%TT = 0.5*(TT + TT');

TT_dtp = kron(R, conj(R)) + kron(Q, eye(D)) + kron(eye(D), conj(Q)) + dx*kron(OP, conj(OP));
TT_dtm = kron(R, conj(R)) + kron(Q, eye(D)) + kron(eye(D), conj(Q)) - dx*kron(OP, conj(OP));
TTeig = sort((eig(TT)));
TTeig_dtp = sort((eig(TT_dtp)));
TTeig_dtm = sort((eig(TT_dtm)));

%dTTeig = 0.5*( (TTeig_dtp - TTeig)/dt - (TTeig - TTeig_dtm)/dt );
dTTeig = (TTeig_dtp - TTeig)/dx;
    
dcorr = sum(dTTeig);

dcorr_std = dTTeig(1);

corr = trace(transpose(matl)*OP*matr*(OP'));



