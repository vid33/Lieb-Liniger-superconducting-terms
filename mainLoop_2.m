matl = eye(D);
l = reshape(eye(D), 1, D^2); 

precision = 1e-14;
        
       % [matr, r_zero_out] = lrEigenvector_fast(R,Q,r_zero, D, precision);
        
        [matr, delta_r_zero] = lrEigenvector_delta(R, Q, matr, delta_r_zero, D, precision);
    
        
        [F, matF ] = calculateF( R, Rkin, Q, matr, potential, uu, interaction, Fzero, D, precision);
        matF_old = matF;


    matr = (1/trace(matr))*matr; %r = reshape(matr, 1, D^2);
        
    calculateEnergy;
    

        deltaEnergy = energyDensity_new - energyDensity_old; 
        

%    Rstar = matr*R'/matr;
 %   Rtangent = transpose(matF)*R-R*transpose(matF)...
 %       +(1/(2*mass))*((Rkin*Q - Q*Rkin)+ ((Rkin*R - R*Rkin)*matr*R')/matr...
 %       + (R*R' - R'*R)*Rkin)...
 %       +potential*R + uu*((matr*R')/matr + R')...
 %       +interaction*((R*R*matr*R')/matr+R'*R*R);
     Rstar = (matr*R')/matr;
    Rtangent = transpose(matF)*R-R*transpose(matF)...
        +(1/(2*mass))*((Rkin*Q - Q*Rkin)+ (Rkin*R - R*Rkin)*Rstar...
        + (R*R' - R'*R)*Rkin)...
        +potential*R + uu*(Rstar + R')...
        +interaction*(R*R*Rstar + R'*R*R);
    
    
    normRtangent = trace(Rtangent*matr*Rtangent');
  
    %*<VW | H | Psi(R,Q) > 
    gradOverlap = trace( transpose(matF)*matr*Qtangent' +transpose(matF)*R*matr*Rtangent'...
        +(1/(2*mass))*Rkin*matr*((Rtangent'*Q' - Q'*Rtangent')+(R'*Qtangent' - Qtangent'*R'))...
        +potential*R*matr*Rtangent' ...
        + uu*matr*(R'*Rtangent' + Rtangent'*R') ...
        + interaction*R*R*matr*(R'*Rtangent'+Rtangent'*R') ); 
    
    deltaEnergyLinear = -1*dt*(gradOverlap+conj(gradOverlap));
        
        info;




