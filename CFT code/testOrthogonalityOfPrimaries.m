clear;

D=8;


Lambda = sqrt(20);
pi_prefactor = (1i/2)*sqrt(2*Lambda);

% for the vertex operator : exp(i beta \Phi) :
beta = 1;

uu = -5;
potential = 10;
interaction = 0;
fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction); 
load(fIn);


dPhi_left = l*(1/sqrt(2*Lambda)*( kron(Rkin, eye(D)) + kron(eye(D), conj(Rkin)) ) ...
                +pi_prefactor*( kron(R, eye(D)) - kron(eye(D), conj(R)) )  );

ddPhi_left = pi_prefactor*l*(kron(Rkin, eye(D)) -  kron(eye(D), conj(Rkin)) );
ddPhi_right = pi_prefactor*(kron(Rkin, eye(D)) -  kron(eye(D), conj(Rkin)) )*r;


H_right = ( (1/Lambda)*kron(Rkin, conj(Rkin)) ...
                                - ((Lambda^2)/4)*( kron(R*R, eye(D)) + kron(eye(D), conj(R*R)) ) ...
                                + ((Lambda^2)/2)*( kron(R, conj(R)) ) )*r;
                            
vertexOp_right = ( expm( (-beta/sqrt(2*Lambda) )*1i*( kron(R, eye(D)) +  kron(eye(D), conj(R))  ) ) )*r;
   


                            
dPhi_H_overlap = dPhi_left*H_right;

dPhi_ddPhi_overlap = dPhi_left*ddPhi_right;

ddPhi_H_overlap = ddPhi_left*H_right;

dPhi_vertexOp_overlap = dPhi_left*vertexOp_right;


%eKinDensity = (1/(2*mass))*trace(Rkin*matr*Rkin');
%ePotDensity = potential*trace(R*matr*R');
%uuDensity = conj(uu)*trace(matr*R'*R')+ uu*trace(R*R*matr);
%eInteractionDensity = interaction*trace(R*R*matr*R'*R');
%energyDensity_new = eKinDensity+  ePotDensity + uuDensity+ eInteractionDensity;
%energyDensity_effective = (eKinDensity+eInteractionDensity)*(1/(occupationDensity^3));
%interaction_effective = interaction/occupationDensity; 


HHmat = (1/(2*mass))*kron(Rkin, conj(Rkin)) + potential*kron(R, conj(R)) ...
    + conj(uu)*kron(eye(D), conj(R*R))+ uu*kron(R*R, eye(D));


testIfEigvec = (dPhi_left*HHmat)./dPhi_left;

