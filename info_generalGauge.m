fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
        fprintf('D = %d, v = %d, uu  = %d, c = %d \n', D, potential, uu, interaction);
        fprintf('AT STEP %d\n', currentStep);
        fprintf('Delta energy to linear order is %d \n' , deltaEnergyLinear );
        fprintf('Delta energy is %d \n' , deltaEnergy );
        fprintf('Occupation density is %d \n' ,occupationDensity );
        fprintf('Ekin density is %d \n' ,eKinDensity );
        fprintf('Epot density is %d \n' ,ePotDensity );
        fprintf('Einteraction density is %d \n' ,eInteractionDensity );
        fprintf('Total energy density is %6.10e \n' ,energyDensity_new );
        fprintf('Energy density ignoring chemical pot. %d\n' , (energyDensity_new - ePotDensity) );
        fprintf('Effective energy density %6.10e at effective c %6.10e\n', energyDensity_effective, interaction_effective);
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');   
 
        fprintf('Delta T is %d \n', dt);
    %    fprintf('The state norm is %d \n', l*r);
        %fprintf('Initial norm is %d \n', initial_norm);
        fprintf('Norm of the gradient is %d \n', normRtangent);
        
        fprintf('delta energy is %d \n', deltaEnergy);
      %  fprintf('delta energy minus delta energy linear %d \n', deltaEnergy-deltaEnergyLinear);
        fprintf('delta energy over delta energy linear %d \n', deltaEnergy/deltaEnergyLinear);
       
        fprintf('min svd(matr) is %d\n', min(svd(matr)));
        fprintf('min svd(matl) is %d\n', min(svd(matl)));
        fprintf('min schmidt is %d \n', min(svd(matl_t^(1/2)*matr^(1/2))) );


        
    
        fprintf('%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%\n');
        
        