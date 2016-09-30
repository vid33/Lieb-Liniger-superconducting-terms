text_info = true;
graphic_info = false;

%for calculating F and density matrices
bicg_maxiter=10000;
bicg_tolerance=1e-12;

if graphic_info == true
    figureName = sprintf('D=%d, v=%d, c=%d', D, potential, interaction);
    energyPlot = figure('Name', figureName); 
    figure(energyPlot);
    set(0,'CurrentFigure', energyPlot);
    xlabel('time');
    ylabel('energy');
end

%Create a stop button so loop below exits nicely
figh = figure;

global IS_ABORTED;
IS_ABORTED = false;

btn = uicontrol('style', 'pushb', 'string', 'Abort', ...
                'callback', @doAbort);
drawnow;
            
while 1  
   
    if (currentStep == 1)  
        deltaEnergy = 0;
    end
    if exist('r_zero', 'var') == 0 
        r_zero= cpxrand(D,D)/sqrt(D);
        r_zero = r_zero+r_zero';
        r_zero = reshape(transpose(r_zero), D^2, 1);
    end
    if VARY_Q == true
        if exist('r_zero_eigs', 'var') == 0 || exist('l_zero_eigs', 'var') == 0 
            r_zero_eigs= cpxrand(D,D)/sqrt(D);
            r_zero_eigs = r_zero_eigs+r_zero_eigs';
            r_zero_eigs = reshape(transpose(r_zero_eigs), D^2, 1);
            l_zero_eigs = reshape(transpose(r_zero_eigs), D^2, 1);
        end
    end
    if exist('F_zero', 'var') == 0 
        F_zero = zeros(D^2, 1);
    end
    %end 
    
    calculateEnergy_generalGauge;
    
    if (currentStep > 1)
        deltaEnergy = energyDensity_new - energyDensity_old;
    end       
    
    Rstar_r = matr*(R'/matr); Rstar_l = (matl_t\R')*matl_t; 
    Qstar_r = matr*(Q'/matr); Qstar_l = (matl_t\Q')*matl_t; 
    
    Rtangent = (matl_t\transpose(matF))*R-R*(matl_t\transpose(matF))...
        +(1/(2*mass))*(Qstar_l*Rkin - Rkin*Qstar_r + R*( Rstar_l*Rkin - Rkin*Rstar_r))...
        +potential*R + conj(uu)*(Rstar_r + Rstar_l) ...
        + interaction*(Rstar_l*(R^2) + (R^2)*Rstar_r);
    
    Ktangent = (1i/2)*(Rtangent'*matl_t*R - R'*matl_t*Rtangent);
    Qtangent = -Rstar_l*Rtangent;
    
    normRtangent = trace(matl_t*Rtangent*matr*Rtangent');
    
    gradOverlap = trace( transpose(matF)*matr*Qtangent' +transpose(matF)*R*matr*Rtangent'...
        +(1/(2*mass))*matl_t*Rkin*matr*( Rtangent'*Q' - Q'*Rtangent' + R'*Qtangent' - Qtangent'*R')...
        +potential*matl_t*R*matr*Rtangent' ...
        + uu*matl_t*matr*(R'*Rtangent' + Rtangent'*R')  ...
        + interaction*matl_t*(R^2)*matr*(R'*Rtangent'+Rtangent'*R') ); 
    
    deltaEnergyLinear = -1*dt*(gradOverlap+conj(gradOverlap));
    
    energyDensity_old = energyDensity_new;
    
    R = R - dt*Rtangent;
    
    if VARY_Q == true
        Q = Q - dt*Qtangent;
    elseif VARY_Q == false
        K = K - dt*Ktangent;
        K = (1/2)*(K + K');
    
        %Rstar_l = matl_t\(R'*matl_t); 
    %    Rstar_l = matl_t_inv*(R'*matl_t);
     % Q = (-1/2)*Rstar_l*R - 1i*(matl_t\K);
     Q = matl_t\(-(1/2)*R'*matl_t*R - 1i*K);
    end
    
    if VARY_Q == false
        [ matr, r_zero ] = calculate_r_fast( R, Q, r_zero, D, bicg_maxiter, bicg_tolerance);
    elseif VARY_Q == true
       [ matl, matr, l_zero_eigs, r_zero_eigs, l_eigval, r_eigval] = calculate_lr_eigs( R, Q, l_zero_eigs, r_zero_eigs, D, bicg_maxiter, bicg_tolerance);
       Q = Q - (r_eigval/2)*eye(D);
       fprintf('r_eigval before changing Q is %d\n', r_eigval);
        fprintf('l_eigval before changing Q is %d\n', l_eigval);
        [ matl, matr, l_zero_eigs, r_zero_eigs, l_eigval, r_eigval] = calculate_lr_eigs( R, Q, l_zero_eigs, r_zero_eigs, D, bicg_maxiter, bicg_tolerance);
       fprintf('r_eigval after changing Q is %d\n', r_eigval);
        fprintf('l_eigval after changing Q is %d\n', l_eigval);
       matl_t = transpose(matl);
       %break;
    end
    tmp_norm = trace(matl_t*matr); matr = matr/tmp_norm;

    Rkin = Q*R - R*Q;
    
    if SYMMETRIC_GAUGE == true
        symmetricGauge;
    end
    
    
    [F] = calculateFL( R, Rkin, Q, matr, matl, potential, uu, interaction, F_zero, D, bicg_maxiter, bicg_tolerance );
    F_zero = transpose(F);
    matF = transpose(reshape(F, D, D));
    
            
    if rem(currentStep,1)==0      
        if text_info == true
           info_generalGauge;
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
        close(figh);
        break;
    end  
    
    drawnow;
    if IS_ABORTED == true
        close(figh);
        break;
    end
    
    currentStep = currentStep+1; currentTime = currentTime + dt;
    
end


