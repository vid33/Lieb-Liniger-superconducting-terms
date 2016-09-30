%(psi | psidag psidag ) overlap  
clear;

D=32; potential= -1; interaction=1;

fIn = sprintf('data/D=%d/cpxD=%dv=%dc=%d.mat', D , D, -potential, interaction);
load(fIn);

p1 = 2;  p2 = 3;

p = p1+p2;

proj = eye(D^2,D^2) - kron(r, l);

Tp1_plus = - kron(Q, eye(D,D)) - kron( eye(D,D), conj(Q)) - kron(R, conj(R)) - kron(r,l) + 1i*p1*eye(D^2, D^2);
Tp1_minus = - kron(Q, eye(D,D)) - kron( eye(D,D), conj(Q)) - kron(R, conj(R)) - kron(r,l) - 1i*p1*eye(D^2, D^2);

Tp2_plus = - kron(Q, eye(D,D)) - kron( eye(D,D), conj(Q)) - kron(R, conj(R)) - kron(r,l) + 1i*p2*eye(D^2, D^2);
Tp2_minus = - kron(Q, eye(D,D)) - kron( eye(D,D), conj(Q)) - kron(R, conj(R)) - kron(r,l) - 1i*p2*eye(D^2, D^2);

Tp_minus = - kron(Q, eye(D,D)) - kron( eye(D,D), conj(Q)) - kron(R, conj(R)) - kron(r,l) - 1i*p*eye(D^2, D^2);
    

Reye = kron(R, eye(D));
eyeR = kron(eye(D), conj(R));
  
overlap  = (l*(Reye*proj)/Tp_minus)*proj*((eyeR*proj)/Tp1_minus)*proj*eyeR*r ...
    +(l*(Reye*proj)/Tp_minus)*proj*((eyeR*proj)/Tp2_minus)*proj*eyeR*r ...
    +(l*(eyeR*proj)/Tp2_plus)*proj*((Reye*proj)/Tp1_minus)*proj*eyeR*r ...
    +(l*(eyeR*proj)/Tp1_plus)*proj*((Reye*proj)/Tp2_minus)*proj*eyeR*r ...
    +(l*(eyeR*proj)/Tp2_minus)*proj*((eyeR*proj)/Tp_minus)*proj*Reye*r ...
    +(l*(eyeR*proj)/Tp1_minus)*proj*((eyeR*proj)/Tp_minus)*proj*Reye*r; 

fprintf('Overlap is %g + i%g\n', real(overlap), imag(overlap));