deltaEnergy = 0;

halt_loop = 0;

wait = 0;
step = 1;


%calculate objects rescale recalculate objects
T = kron(Q, eye(D))+kron(eye(D), conj(Q))+kron(R, conj(R));
[l, r] = calculatelr_slow(T, D);
[matl, matr] = calculateMatlr(l, r, D);

Rkin = commutator(Q, R);
    
occupationDensity = (1/(l*r))*l*kron(R, conj(R))*r;
eKinDensity = (1/(l*r))*(1/(2*mass))*l*kron(Rkin, conj(Rkin))*r;
ePotDensity = (1/(l*r))*potential*l*kron(R, conj(R))*r;
eInteractionDensity = (1/(l*r))*interaction*l*kron(R*R, conj(R)*conj(R))*r;
energyDensity_old = eKinDensity+  ePotDensity + eInteractionDensity;
energyDensity_effective = (eKinDensity+eInteractionDensity)*(1/(occupationDensity^3));

interaction_effective = interaction/occupationDensity;  

F = (-1*l*((1/(2*mass))*kron(commutator(Q,R), commutator(conj(Q),conj(R))) ...
    + potential*kron(R, conj(R))+interaction*kron(R^2, conj(R)^2))    ...
    +energyDensity_old*l )/(T + kron(r, l ));
matF = transpose(reshape(F, D, D));

Rtangent = commutator(transpose(matF),R)+(1/(2*mass))*(commutator(Rkin, Q)+ commutator(Rkin, R)*(matr*(hconj(R))*(matr^(-1)))...
        + commutator(R, hconj(R))*Rkin)...
        +potential*R + interaction*(R*R*matr*hconj(R)*(matr^(-1))+hconj(R)*R*R);
 
Ktangent = i*0.5*(hconj(Rtangent)*R-hconj(R)*Rtangent);



while halt_loop  == 0   
    
    Rkin = commutator(Q, R); 
    
    Qtangent = -hconj(R)*Rtangent;
    
     gradOverlap = F*(kron(eye(D,D), conj(Qtangent))+ kron(R, conj(Rtangent)))*r...
    + (1\(2*mass))*l*(kron(Rkin, (commutator(conj(Q), conj(Rtangent)) +commutator(conj(Qtangent), conj(R)))  ))*r...
    + potential*l*kron(R, conj(Rtangent))*r + interaction*l*kron(R^2, conj(Rtangent)*conj(R))*r ...
    + interaction*l*kron(R^2, conj(R)*conj(Rtangent))*r ;
    
    deltaEnergyLinear = -1*dt*(1/(l*r))*(gradOverlap+conj(gradOverlap));

    if step > 10
        normRtangent_old = normRtangent;
    end
    
    
    
    normRtangent = (1/(l*r))*l*(kron(Rtangent, conj(Rtangent)))*r;
    
    
    K = K-dt*Ktangent;
    R = R - dt*Rtangent;      
    Q = -0.5*hconj(R)*R-i*K;
    
    %rescale to get desnity 1
   % rho = (1/(l*r))*l*kron(R, conj(R))*r;
   % R = sqrt((1/rho))*R; Q = (1/rho)*Q; K = (1/rho)*K;

    
    Rkin = commutator(Q, R);
    T = kron(Q, eye(D))+kron(eye(D), conj(Q))+kron(R, conj(R));
    [l, r] = calculatelr_slow(T, D);
    [matl, matr] = calculateMatlr(l, r, D);
    
    occupationDensity = (1/(l*r))*l*kron(R, conj(R))*r;
    eKinDensity = (1/(l*r))*(1/(2*mass))*l*kron(Rkin, conj(Rkin))*r;
    ePotDensity = (1/(l*r))*potential*l*kron(R, conj(R))*r;
    eInteractionDensity = (1/(l*r))*interaction*l*kron(R*R, conj(R)*conj(R))*r;
    energyDensity_new = eKinDensity+  ePotDensity + eInteractionDensity;
    energyDensity_effective = (eKinDensity+eInteractionDensity)*(1/(occupationDensity^3));
    
    if(step > 1)
        deltaEnergy = energyDensity_new - energyDensity_old;
    end
    
    interaction_effective = interaction/occupationDensity;  
    
    
    F = (-1*l*((1/(2*mass))*kron(commutator(Q,R), commutator(conj(Q),conj(R))) ...
    + potential*kron(R, conj(R))+interaction*kron(R^2, conj(R)^2))    ...
    +0*energyDensity_new*l )/(T + kron(r, l) );
    %F = -1*energyDensity_new*l;
    matF = transpose(reshape(F, D, D));
    


    Rtangent = commutator(transpose(matF),R)+(1/(2*mass))*(commutator(Rkin, Q)+ commutator(Rkin, R)*(matr*(hconj(R))*(matr^(-1)))...
        + commutator(R, hconj(R))*Rkin)...
        +potential*R + interaction*(R*R*matr*hconj(R)*(matr^(-1))+hconj(R)*R*R);
    
    Ktangent = i*0.5*(hconj(Rtangent)*R-hconj(R)*Rtangent);

    
    if rem(step,50)==0
        
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
        fprintf('AT STEP %d\n', step);
        fprintf('Delta energy to linear order is %d \n' , deltaEnergyLinear );
        fprintf('Delta energy is %d \n' , deltaEnergy );
        fprintf('Occupation density is %d \n' ,occupationDensity );
        fprintf('Ekin density is %d \n' ,eKinDensity );
        fprintf('Epot density is %d \n' ,ePotDensity );
        fprintf('Einteraction density is %d \n' ,eInteractionDensity );
        fprintf('Total energy density is %d \n' ,energyDensity_new );
        fprintf('Energy density ignoring chemical pot. %d\n' , (energyDensity_new - ePotDensity) );
        fprintf('Effective energy density %d at effective c %d\n', energyDensity_effective, interaction_effective);
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');   
 
        fprintf('Delta T is %d \n', dt);
        fprintf('The state norm is %d \n', l*r);
        fprintf('Norm of the gradient is %d \n', normRtangent);
        fprintf('delta energy over delta energy linear %d \n', deltaEnergy/deltaEnergyLinear);
        
    
        fprintf('%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%\n');
        
        
    end
    
    if normRtangent< tolerance
        fprintf('Exiting and saving....\n');
        fprintf('Effective energy density is %d, at effective c %d\n',  energyDensity_effective, interaction_effective);
        %save(fOut);
        break;
    end
    
    if (step>500 && (deltaEnergy/deltaEnergyLinear < 1.05 || deltaEnergy/deltaEnergyLinear > 0.95) && normRtangent < 7e-6)
       % dt = 0.0075;
    end
    
     if (step>500 && (deltaEnergy/deltaEnergyLinear > 1.2 || deltaEnergy/deltaEnergyLinear < 0.8))
    
        fprintf('Warning something is ffing wrong!....\n');
        save('data/bad_exit.mat');
        %break;
    end
    
%     if (step>200 && (deltaEnergy/deltaEnergyLinear > 1.1 || deltaEnergy/deltaEnergyLinear < 0.9))
%         if(dt > tolerance)
%             dt = dt*0.5;
%         else
%             fprintf('Exiting and saving....\n');
%             fprintf('Effective energy density is %d, at effective c %d\n',  energyDensity_effective, interaction_effective);
%            % save(fOut);
%             break;
%         end
%      end
%     if (step>200 && (deltaEnergy/deltaEnergyLinear < 1.05 || deltaEnergy/deltaEnergyLinear > 0.95))
%         if(dt < 0.01)
%             dt = dt/0.9;
%         end
%     end
    

    
    
    energyDensity_old = energyDensity_new;
    step = step+1;
    
end


