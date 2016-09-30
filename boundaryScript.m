
options=odeset('RelTol',1e-6, 'AbsTol', 1e-9);

fOut = sprintf('data_LL/D=%d/boundary/stripD=%dL=%dN=%duu=%dv=%dc=%d.mat', D, D, L, N, uu, potential, interaction);

TDVP = false;

haltITE = 0;

currentStep = 0;

leftStep = true;

while haltITE==0
    currentStep = currentStep+1;
    fprintf('Current step: %d\n', currentStep);

    %T_action_right_const(t_in, x_in, R, Q, D)
    solution_l = ode45(@(t,x) T_action_right_const(t, x, R, Q, D), [0 L], transpose(vl), options);

    position_l_ode = solution_l.x;
    l_ode = solution_l.y;

    l = deval(solution_l,position);

    for kk=1:N
        matl_cell{kk} = transpose(reshape(transpose(l(:, kk)), D, D));
    end
    
    norm_tmp = trace(transpose(matl_cell{end})*matvr);
    matvr = matvr/norm_tmp;
    vr = vr/norm_tmp;
   
    %function out=T_action_left_const(t_in, x_in, R, Q, D)
    solution_r = ode45(@(t,x) T_action_left_const(t, x, R, Q, D), [L 0], vr, options);

    position_r_ode = solution_r.x;
    r_ode = solution_r.y;
    
    fprintf('Numel position_r_ode is %d\n', numel(position_r_ode));

    r = deval(solution_r, position);

    for kk=1:N
            matr_cell{kk} = transpose(reshape((r(:, kk)), D, D));
            testNorm(kk) = trace(transpose(matl_cell{kk})*matr_cell{kk});
    end
    
    


    %F_action_right_const(t_in, F_in, R, Q, Rkin, mass, potential, uu, interaction, D, solution_l)
    solution_FL = ode45(@(t,x) F_action_right_const(t, x, R, Q, Rkin, mass, potential, uu, interaction, D, solution_l),...
        [0 L], zeros(1, D^2), options);
    
    FL = deval(solution_FL,position);

    for kk=1:N
        matFL{kk} = transpose(reshape((FL(:, kk)), D, D)); 
    end
    
    %F_action_left_const(t_in, F_in, R, Q, Rkin, mass, potential, uu, interaction, D, solution_r)
    solution_FR = ode45(@(t,x) F_action_left_const(t, x, R, Q, Rkin, mass, potential, uu, interaction, D, solution_r),...
        [L 0], zeros(D^2, 1), options);
    
    FR = deval(solution_FR, position);

    for kk=1:N
        matFR{kk} = transpose(reshape(transpose(FR(:, kk)), D, D)); 
    end
    
    for kk=1:N
        energyDensityL(kk) = trace(transpose(matFL{kk})*matr_cell{kk});
        energyDensityR(kk) = trace(transpose(matl_cell{kk})*matFR{kk});
        
        matl_t = transpose(matl_cell{kk});
        
        occupationDensity(kk) = trace(matl_t*R*matr_cell{kk}*R');
        eKinDensity(kk) = trace(matl_t*Rkin*matr_cell{kk}*Rkin');
        eInteractionDensity(kk) = trace(matl_t*R*R*matr_cell{kk}*R'*R');
        
        energyDensity(kk) = (1/(2*mass))*eKinDensity(kk) ...
            + potential*occupationDensity(kk)  ...
            + uu*trace(matl_t*matr_cell{kk}*R'*R')+ uu*trace(matl_t*R*R*matr_cell{kk}) ...
            + interaction*eInteractionDensity(kk);
        
        Uschmidt = transpose(matl_cell{kk}^0.5)*(matr_cell{kk}^0.5);
        schmidt_sq(:, kk) = svd(Uschmidt*Uschmidt');
        
        entropy(kk) = sum(-schmidt_sq(:, kk).*log(schmidt_sq(:, kk)));
    end
    energyTotal = trapz(position, energyDensity);
    
    
    fprintf('Total energy is %d\n', energyTotal);
    
   if rem(currentStep,250)==0
        save('data_LLuu/tmp_save_boundary.mat');
   end 
   
    
    WR = inv(matl_cell{N})*matFL{N}*VRsqrt;
    WR = WR - ((VRsqrt')*WR)*VRsqrt;
    if TDVP==true
        VRsqrt = VRsqrt - dt*WR;
    else
        if leftStep == true 
            [VmatFL, DmatFL] = eig(inv(matl_cell{N})*matFL{N});
            [DmatFL, DmatFL_index] = sort(diag(DmatFL));
            VRsqrt = VmatFL(:, DmatFL_index(1));
        end
    end
    matvr  = kron(VRsqrt', VRsqrt);
    matvr = conj(matvr);
    vr = reshape(transpose(matvr), D^2, 1);

    WL = VLsqrt*transpose(matFR{1})*transpose(inv(matr_cell{1}));
    WL = WL - (WL*(VLsqrt'))*VLsqrt;
    
    if TDVP==true
        VLsqrt  = VLsqrt - dt*WL;
    else
        if leftStep == false
            [VmatFR, DmatFR] = eig(inv(matr_cell{1})*matFR{1});
            [DmatFR, DmatFR_index] = sort(diag(DmatFR));
            VLsqrt = VmatFR(:, DmatFR_index(1));
            VLsqrt  = transpose(VLsqrt);
        end
    end
    matvl  = kron(VLsqrt', VLsqrt);
    vl = reshape(transpose(matvl), 1, D^2);
    %norm = trace(transpose(matvl)*matr);
    %matvl = matvl/norm;
    %vl = vl/norm;
    %vlsqrt = vlsqrt/sqrt(norm);
   
    leftStep = mod(leftStep +1, 2);   
    
    WRmat = kron(WR', WR);
    WRmat = conj(WRmat);
    WLmat = kron(WL', WL);
    norm_WR = trace(transpose(matl_cell{end})*WRmat);
    norm_WL = trace(transpose(WLmat)*matr_cell{1});
    
    fprintf('Norm of WL %d, Norm of WR %d\n', norm_WL, norm_WR);
       
end
    
    
