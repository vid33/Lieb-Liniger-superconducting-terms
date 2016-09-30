occupationDensity = trace(R*matr*R');
eKinDensity = (1/(2*mass))*trace(Rkin*matr*Rkin');
ePotDensity = potential*trace(R*matr*R');
uuDensity = conj(uu)*trace(matr*R'*R')+ uu*trace(R*R*matr);
eInteractionDensity = interaction*trace(R*R*matr*R'*R');
energyDensity_new = eKinDensity+  ePotDensity + uuDensity+ eInteractionDensity;
energyDensity_effective = (eKinDensity+eInteractionDensity)*(1/(occupationDensity^3));
interaction_effective = interaction/occupationDensity; 


