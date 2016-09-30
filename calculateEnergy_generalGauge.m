occupationDensity = trace(matl_t*R*matr*R');
eKinDensity = (1/(2*mass))*trace(matl_t*Rkin*matr*Rkin');
ePotDensity = potential*trace(matl_t*R*matr*R');
uuDensity = conj(uu)*trace(matl_t*matr*R'*R')+ uu*trace(matl_t*R*R*matr);
eInteractionDensity = interaction*trace(matl_t*R*R*matr*R'*R');
energyDensity_new = eKinDensity+  ePotDensity + uuDensity+ eInteractionDensity;
energyDensity_effective = (eKinDensity+eInteractionDensity)*(1/(occupationDensity^3));
interaction_effective = interaction/occupationDensity; 


