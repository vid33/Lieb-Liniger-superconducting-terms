text_info = true;
graphic_info = false;

%for calculating F and density matrices
bicg_maxiter=20000;
bicg_tolerance=1e-14;

if graphic_info == true
    figureName = sprintf('D=%d, v=%d, c=%d', D, potential, interaction);
    energyPlot = figure('Name', figureName); 
    figure(energyPlot);
    set(0,'CurrentFigure', energyPlot);
    xlabel('time');
    ylabel('energy');
end


while 1  
   
    if (currentStep == 1)  
        deltaEnergy = 0;   
        r_zero= cpxrand(D,D)/sqrt(D);
        r_zero = r_zero+r_zero';
        r_zero = reshape(transpose(r_zero), D^2, 1);
        F_zero = zeros(D^2, 1);
    end
    
    if D ~=1
  %      [l, r, matr, initial_norm, r_zero_out] = lrEigenvector_fast(R,Q,r_zero, D);
        [ matr, r_zero ] = calculate_r_fast( R, Q, r_zero, D, bicg_maxiter, bicg_tolerance);
       
    end
    
    calculateEnergy;
    
    if (currentStep > 1)
        deltaEnergy = energyDensity_new - energyDensity_old;
    end       
        
    %calculateFL( R, Rkin, Q, matr, matl, potential, uu, interaction, F_zero, D, maxiter, tolerance )
    [F] = calculateFL( R, Rkin, Q, matr, matl, potential, uu, interaction, F_zero, D, bicg_maxiter, bicg_tolerance );
    F_zero = transpose(F);
    matF = transpose(reshape(F, D, D));
    
    Rstar = matr*R'/matr;
%    Rtangent = transpose(matF)*R-R*transpose(matF)...
%        +(1/(2*mass))*((Rkin*Q - Q*Rkin)+ (Rkin*R - R*Rkin)*Rstar...
%        + (R*R' - R'*R)*Rkin)...
%        +potential*R + interaction*(R*R*Rstar+R'*R*R);
 
    Rtangent = transpose(matF)*R-R*transpose(matF)...
        +(1/(2*mass))*((Rkin*Q - Q*Rkin)+ (Rkin*R - R*Rkin)*Rstar...
        + (R*R' - R'*R)*Rkin)...
        +potential*R + conj(uu)*(Rstar + R')...
        +interaction*(R*R*Rstar + R'*R*R);
    
    Ktangent = 1i*0.5*(Rtangent'*R-R'*Rtangent);
    Qtangent = -R'*Rtangent;
   
    normRtangent = trace(Rtangent*matr*Rtangent');
    
    gradOverlap = trace( transpose(matF)*matr*Qtangent' +transpose(matF)*R*matr*Rtangent'...
        +(1/(2*mass))*Rkin*matr*((Rtangent'*Q' - Q'*Rtangent')+(R'*Qtangent' - Qtangent'*R'))...
        +potential*R*matr*Rtangent' ...
        + conj(uu)*matr*(R'*Rtangent' + Rtangent'*R') ...
        + interaction*R*R*matr*(R'*Rtangent'+Rtangent'*R') ); 
    
    deltaEnergyLinear = -1*dt*(gradOverlap+conj(gradOverlap));
    
        
    R = R - dt*Rtangent;  
    K = K - dt*Ktangent;
    K = 0.5*(K + K');
    
    Q = -0.5*(R')*R-1i*K;
    Rkin = Q*R - R*Q;
    
    energyDensity_old = energyDensity_new;
    currentStep = currentStep+1; currentTime = currentTime + dt;
    
    if rem(currentStep,1)==0      
        if text_info == true
           info;
        end
    end
    if graphic_info == true
        figure(energyPlot);
        if rem(currentStep,100) == 0;
            clf(energyPlot);
            xlabel('time');
            ylabel('energy');
        end
        hold on;
        plot(currentTime, real(energyDensity_new), 'x');
        drawnow;
        hold off;
    end
    
    if rem(currentStep,250)==0
        save('data_LLuu/tmp_save.mat');
    end  
        
    if normRtangent< tolerance
        fprintf('Exiting and saving....\n');
        fprintf('Effective energy density is %d, at effective c %d\n',  energyDensity_effective, interaction_effective);
       %save(fOut);
        break;
    end  
end


