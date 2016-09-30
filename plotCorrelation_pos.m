clear;

D=16;

uu = -5;
potential = 10;
interaction = 0;

interpolate = false;
useExpm = false;  %Do the calculation the long way, D^6, to check stuff
plotDominant = false; %Plot the part of correlator coming just from the 

field = false;
PI = false; %not pi, to avoid conflict with 3.1415..
dPI = false;
dPhi = false;
hamiltonian = false;
density = false;
ekin = false;
ekin_density = false;
dPhi_hamiltonian = false;
vertexOp = false;
dPhi_vertexOp = false;
vertexOp_hamiltonian = false;
dPhi_dPI = true;

Lambda = sqrt(20);
pi_prefactor = (1i/2)*sqrt(2*Lambda);

% for the vertex operator : exp(i beta \Phi) :
beta = 0.1;

fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction); 
load(fIn);
Rkin = Q*R - R*Q;

%load correlation length values
fIn= sprintf('data_LLuu/lr_eigenvectors/uu=%dv=%dc=%deigenvector=%d.mat',uu,  potential, interaction , 2);
load(fIn, 'mu', 'dim');

%set maximum x value in terms of correlation length
scale = 4;
for kk=1:numel(dim)
    if (dim(kk) == D)
        x_max = real(-1/mu(kk))*scale;
    end
end

%For hamiltonian - hamiltonian correlator
eDensity_rescaled = (1/Lambda)*trace(Rkin*matr*Rkin') ...
      - ((Lambda^2)/4)*( trace(R*R*matr) + trace(matr*R'*R') ) ...
          + ((Lambda^2)/2)*trace(R*matr*R') ;

x_min = 0;

options=odeset('RelTol',1e-10, 'AbsTol', 1e-13);

if (field ==true)
    [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), [x_min x_max], reshape(R, D^2, 1), options );
elseif (PI ==true)
    [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), [x_min x_max], ...
                                                    reshape(pi_prefactor*(R - R'), D^2, 1), options );
elseif (dPI ==true )
    [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), [x_min x_max], ...
                                                    reshape(pi_prefactor*(Rkin - Rkin'), D^2, 1), options );
elseif (dPhi ==true || dPhi_hamiltonian == true || dPhi_vertexOp == true || dPhi_dPI==true )
    [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), [x_min x_max], ...
                                       reshape( 0.5*(1/sqrt(2*Lambda)*(Rkin + Rkin') + pi_prefactor*(R - R') ), D^2, 1), options );
elseif (hamiltonian == true )
    [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), [x_min x_max], ...
                reshape( (1/Lambda)*(Rkin')*Rkin ...
                -((Lambda^2)/4)*(R*R + R'*R') + ((Lambda^2)/2)*R'*R, D^2, 1 ) );
elseif (density == true)
    [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), [x_min x_max], reshape(R'*R, D^2, 1), options );
elseif (ekin == true)
    [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), [x_min x_max], reshape(Rkin, D^2, 1), options );
elseif (ekin_density == true)
    [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), [x_min x_max], reshape(Rkin'*Rkin, D^2, 1), options );
elseif (vertexOp == true || vertexOp_hamiltonian == true)
        %compute (exp(beta*R\times 1...) acting on left and right
        %eigenvectors of T, and exp..T acting on the l*exp(beta*R\times1...) 
                    
        [pos_out1, out_corr1] = ode113(@(t,x) vertex_action(t, x, beta, Lambda, R, D), [0 1], ...
                reshape(eye(D), D^2, 1), options ) ;
                                
        matrix_left_vertex = transpose(reshape(transpose(out_corr1(end, :)), D, D));
        
        [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), [x_min x_max] , ...
                                    reshape(transpose(matrix_left_vertex), D^2, 1), options);
end

correlationFn = zeros(1, numel(pos_out));

for kk=1:numel(pos_out)
    matrix_left = transpose(reshape(transpose(out_corr(kk, :)), D, D));
    
    if (field == true)
        correlationFn(kk) = trace(matrix_left*conj(R)*transpose(matr)) - trace(R*matr)*trace(matr*R');
    elseif (PI == true)
        correlationFn(kk) = trace(matrix_left*pi_prefactor*( transpose(matr)*transpose(R) - conj(R)*transpose(matr) ));
    elseif (dPI == true || dPhi_dPI == true)
        correlationFn(kk) = trace(matrix_left*pi_prefactor*( transpose(matr)*transpose(Rkin) ...
                        - conj(Rkin)*transpose(matr) ));
    elseif (dPhi == true)
        correlationFn(kk) = trace(matrix_left*0.5*( 1/sqrt(2*Lambda)*( transpose(matr)*transpose(Rkin) + conj(Rkin)*transpose(matr) ) ...
           + pi_prefactor*( transpose(matr)*transpose(R) - conj(R)*transpose(matr) ) ) );
    elseif (hamiltonian == true || dPhi_hamiltonian || vertexOp_hamiltonian == true)
        correlationFn(kk) = trace(matrix_left*( (1/Lambda)*conj(Rkin)*transpose(matr)*transpose(Rkin) ...
                               -((Lambda^2)/4)*( conj(R*R)*transpose(matr) + transpose(matr)*transpose(R*R) ) ...
                               +((Lambda^2)/2)*conj(R)*transpose(matr)*transpose(R) ) ) ...
                               - 0*eDensity_rescaled^2 - trace(matrix_left*transpose(matr))*eDensity_rescaled;
    elseif (density == true)
        correlationFn(kk) = ...
            (1/(occupationDensity^2))*(trace(matrix_left*conj(R)*transpose(matr)*transpose(R)) - occupationDensity^2);
    elseif (ekin == true)
        correlationFn(kk) = trace(matrix_left*conj(Rkin)*transpose(matr)) - trace(Rkin*matr)- trace(matr*Rkin');
    elseif (ekin_density == true)
         correlationFn(kk) = ...
            trace(matrix_left*conj(Rkin)*transpose(matr)*transpose(Rkin)) - trace(Rkin*matr*Rkin')^2;
    elseif (vertexOp == true || dPhi_vertexOp == true || hamiltonian_vertexOp == true)
        %compute (exp(beta*R\times 1...) acting on left and right
        %eigenvectors of T, and exp..T acting on the l*exp(beta*R\times1...) 
        
        [pos_out_right, out_corr_right] = ode113(@(t,x) vertex_action_right(t, x, -beta, Lambda, R, D), [0 1], ...
                reshape(transpose(matr), D^2, 1), options ) ;
        matrix_right_vertex = transpose(reshape(transpose(out_corr_right(end, :)), D, D));
        
        correlationFn(kk)  = trace(matrix_left*transpose(matrix_right_vertex));
            
    end
        
end


correlationFn = transpose(correlationFn);

figureName = sprintf('D=%d, c_{eff}=%d', D, interaction_effective);

%figure('Name', figureName); plot(plot_x, plot_y, 'x'); 
%title(figureName);
%xlabel('$x$', 'Interpreter', 'LaTex'); 
%ylabel('$\langle \Psi^\dagger(0) \Psi(x) \rangle$', 'Interpreter', 'LaTex');


figure('Name', figureName); 
%plot( pos_out, sqrt(correlationFn.*conj(correlationFn)) );
%plot( pos_out, (1/2)*(correlationFn + conj(correlationFn)) );

plot(log(pos_out), log( abs(real(correlationFn)) ));

break;

%draw line at correlation length
line([log(x_max/scale); log(x_max/scale)],[ min(log((real(correlationFn))) ); max(log((real(correlationFn))) )], 'Color', 'b');
%draw lines at fractions of correlation length
factor_tmp = 0.75;
line([log(factor_tmp*x_max/scale); log(factor_tmp*x_max/scale)], ...
                [ min(log((real(correlationFn))) ); max(log((real(correlationFn))) )], 'Color', 'g');
factor_tmp = 0.5;
line([log(factor_tmp*x_max/scale); log(factor_tmp*x_max/scale)], ...
                [ min(log((real(correlationFn))) ); max(log((real(correlationFn))) )], 'Color', 'r');
factor_tmp =0.25;
line([log(factor_tmp*x_max/scale); log(factor_tmp*x_max/scale)], ...
                [ min(log((real(correlationFn))) ); max(log((real(correlationFn))) )], 'Color', 'y');

%plot(pos_out, real(correlationFn));
%title(figureName);
xlabel('$\log (x)$', 'Interpreter', 'LaTex'); 
ylabel('$\log ( \langle \Psi^\dagger(0) \Psi(x)  \rangle ) $', 'Interpreter', 'LaTex');

if (useExpm == true)
    
    %delta_x=0.1;
    %plot_x = x_min:delta_x:x_max;
    plot_x = linspace(x_min, x_max, 50);
    
    correlation_x = zeros(1, numel(plot_x));
    %plot_y = zeros(1, numel(plot_x));
    T = kron(Q, eye(D))+ kron( eye(D), conj(Q)) + kron(R, conj(R));

    for k = 1: numel(plot_x)
       
        fprintf('%d ', k);
        if (rem(k,10)==0)
            fprintf('\n');
        end
        
        if (field == true)
            correlation_x(k) = l*(kron(R, eye(D))*expm(plot_x(k)*T)*kron(eye(D), conj(R)))*r;
            %                 - l*kron(R, eye(D))*r*l*kron(eye(D), conj(R))*r;

                %correlation_x(k) = (1/occupationDensity)*l*(kron(R, eye(D,D))*expm(plot_x(k)*T)*kron(eye(D,D), conj(R)))*r;
            %                    - (1/occupationDensity)*l*kron(R, eye(D))*r*l*kron(eye(D), conj(R))*r; 
        
        elseif (PI == true)
            correlation_x(k) = pi_prefactor^2*l*( kron(R, eye(D)) - kron(eye(D), conj(R)) )*expm(plot_x(k)*T)*( ...
                                   kron(R, eye(D)) -kron(eye(D), conj(R)) )*r;
        elseif (dPhi == true)
            correlation_x(k) = (1/4)*l*(1/sqrt(2*Lambda)*( kron(Rkin, eye(D)) + kron(eye(D), conj(Rkin)) ) ...
                +pi_prefactor*( kron(R, eye(D)) - kron(eye(D), conj(R)) ) )*expm(plot_x(k)*T)*( ...
                            1/sqrt(2*Lambda)*( kron(Rkin, eye(D)) + kron(eye(D), conj(Rkin)) ) ...
                            +pi_prefactor*( kron(R, eye(D)) -kron(eye(D), conj(R)) ) )*r;
        elseif (hamiltonian == true)
            correlation_x(k) = l*(  (1/Lambda)*kron(Rkin, conj(Rkin)) ...
                                - ((Lambda^2)/4)*( kron(R*R, eye(D)) + kron(eye(D), conj(R*R)) ) ...
                                + ((Lambda^2)/2)*( kron(R, conj(R)) ) )*expm(plot_x(k)*T)*( ...
                                  (1/Lambda)*kron(Rkin, conj(Rkin)) ...
                                - ((Lambda^2)/4)*( kron(R*R, eye(D)) + kron(eye(D), conj(R*R)) ) ...
                                + ((Lambda^2)/2)*( kron(R, conj(R)) ) )*r - eDensity_rescaled^2;
            
                       % aal_T_disconn(k) = (1/Lambda)*( trace(Rkin*matr*Rkin') ...
       %     - ((Lambda^2)/4)*( trace(R*R*matr) + trace(matr*R'*R') ) ...
       %     + ((Lambda^2)/2)*trace(R*matr*R') ); 
            
        elseif (density == true)
         correlation_x(k) = (1/(occupationDensity^2))*(l*(kron(R, conj(R))*expm(plot_x(k)*T)*kron(R, conj(R)))*r ...
             - occupationDensity^2);
   
        elseif (ekin == true)
         
         correlation_x(k) = l*(kron(Rkin, eye(D,D))*expm(plot_x(k)*T)*kron(eye(D,D), conj(Rkin)))*r ...
                          - trace(Rkin*matr)*trace(matr*Rkin');  
                      
        elseif (ekin_density == true)
        
          correlation_x(k) = l*(kron(Rkin, conj(Rkin))*expm(plot_x(k)*T)*kron(Rkin, conj(Rkin)))*r...
             - l*kron(Rkin, conj(Rkin))*r*l*kron(Rkin, conj(Rkin))*r;
    
        end
    end

    hold all;
    plot(log(plot_x), log(real(correlation_x)), 'x');
    %plot(plot_x, correlation_x, 'x');

    hold off;
end


if (plotDominant == true)
    
    fIn = sprintf('data_LLuu/lr_eigenvectors/uu=%dv=%dc=%deigenvector=%d.mat',uu,  potential, interaction , 2);
        load(fIn, 'dim');
    
    for pp=1:numel(dim)
        if (dim(pp) == D)
            dim_index = pp;
        end
    end
    
    eigenvectors = [ 2 3 4 ]; %eigenvectors to be used
    
    eigenvectorNo_max = 50;
    fIn = sprintf('data_LLuu/lr_eigenvectors/%dEIGVECSuu=%dv=%dc=%d.mat',eigenvectorNo_max, uu,  potential, interaction);
    load(fIn, 'dim', 'mu', 'matl2', 'matr2');
    
    plot_x = linspace(x_min, x_max, 200);
    
    correlationDominant = zeros(1, numel(plot_x));
    
    if (field == true)
        for pp=1: numel(correlationDominant)
            correlationDominant(pp) = trace(matr*R')*trace(R*matr);
        end
    elseif (dPhi == true)
        %for pp=1: numel(correlationDominant)
        %    tmp1 = 0.5*( 1/sqrt(2*Lambda)*( trace(Rkin*matr) + trace(matr*Rkin' ) )  ...
        %            + pi_prefactor*( trace(R*matr) - trace(matr*R') ) );
        %    tmp2 = 0.5*( 1/sqrt(2*Lambda)*( trace(matr*Rkin) + trace(matr*Rkin') ) ...
        %            + pi_prefactor*( trace(matr*R) - trace(matr*R') ) );
        %    correlationDominant(pp) = tmp1*tmp2;
        %end
        %correlationDominant_test = correlationDominant;
    end
       
    
    for pp=1:numel(eigenvectors)
        
        %fIn = sprintf('data_LLuu/lr_eigenvectors/uu=%dv=%dc=%deigenvector=%d.mat',uu,  potential, interaction , eigenvectors(pp));
        %load(fIn, 'mu', 'dim', 'matr2', 'matl2', 'r2', 'l2');
        if (density == true)
            correlationDominant = correlationDominant+ exp(plot_x*mu(dim_index))*trace(matr2{dim_index}*R')*trace(matr2{dim_index}*R')*trace(matr*transpose(matl2{k})*R);
        elseif (PI == true)
               % tmp1 =  pi_prefactor*( trace(R*matr2{dim_index}) -trace(matr2{dim_index}*R') ); 
               % tmp2 = pi_prefactor*( trace(matr*transpose(matl2{dim_index})*R) - trace(matr*R'*transpose(matl2{dim_index})) );
               % correlationDominant = correlationDominant...
               %     + exp(plot_x*mu(dim_index))*tmp1*tmp2;
                
                tmp1 =  pi_prefactor*(  trace(R*matr2{eigenvectors(pp), D-1}) -trace(matr2{eigenvectors(pp), D-1}*R') ); 
                tmp2 = pi_prefactor*(  trace(matr*transpose(matl2{eigenvectors(pp), D-1})*R) ...
                    - trace(matr*R'*transpose(matl2{eigenvectors(pp), D-1}))  );
               correlationDominant =  correlationDominant ...
                    + exp(plot_x*mu{eigenvectors(pp), D-1})*tmp1*tmp2;
        elseif (dPhi == true)
                %tmp1 = 0.5*( 1/sqrt(2*Lambda)*( trace(Rkin*matr2{dim_index}) + trace(matr2{dim_index}*Rkin' ) )  ...
                %    + pi_prefactor*( trace(R*matr2{dim_index}) - trace(matr2{dim_index}*R') ) );
                %tmp2 = 0.5*( 1/sqrt(2*Lambda)*( trace(matr*transpose(matl2{dim_index})*Rkin) + trace(matr*Rkin'*transpose(matl2{dim_index})) ) ...
                %    + pi_prefactor*( trace(matr*transpose(matl2{dim_index})*R) - trace(matr*R'*transpose(matl2{dim_index})) ) );
                %correlationDominant = correlationDominant + exp(plot_x*mu(dim_index) )*tmp1*tmp2;
                
                tmp1 = 0.5*( 1/sqrt(2*Lambda)*( trace(Rkin*matr2{eigenvectors(pp), D-1}) + trace(matr2{eigenvectors(pp), D-1}*Rkin' ) )  ...
                    + pi_prefactor*( trace(R*matr2{eigenvectors(pp), D-1}) - trace(matr2{eigenvectors(pp), D-1}*R') ) );
                tmp2 = 0.5*( 1/sqrt(2*Lambda)*( trace(matr*transpose(matl2{eigenvectors(pp), D-1})*Rkin) + trace(matr*Rkin'*transpose(matl2{eigenvectors(pp), D-1})) ) ...
                    + pi_prefactor*( trace(matr*transpose(matl2{eigenvectors(pp), D-1})*R) - trace(matr*R'*transpose(matl2{eigenvectors(pp), D-1})) ) );
    
                correlationDominant = correlationDominant + exp(plot_x*mu{eigenvectors(pp), D-1})*tmp1*tmp2;
            
        end
    end
    
    clearvars tmp1 tmp2;
    
    hold all;
    plot(log(plot_x), log(real(correlationDominant)));
    %plot(plot_x, correlation_x, 'x');

    hold off;
    
end


if (interpolate == true)
    
    interpolationMax = 0.01:0.005:1.35;
    rangeScale = 0.01;
    interpolationRange = rangeScale*x_max/scale;
    interpolationMidpoint = zeros(1, numel(interpolationMax));
    
    slope = zeros(1, numel(interpolationMax)); intercept = zeros(1, numel(interpolationMax));
    sse = zeros(1, numel(interpolationMax)); rmse = zeros(1, numel(interpolationMax));
    rsquare = zeros(1, numel(interpolationMax)); adjrsquare = zeros(1, numel(interpolationMax));
    confint95 = cell(1, numel(interpolationMax)); 
    confint99 = cell(1, numel(interpolationMax));
    
    fprintf('Performing interpolations \n');
    for zz=1:numel(interpolationMax)
        fprintf('%d ', zz);
        if (rem(zz,10)==0)
           fprintf('\n');
        end
        
        log_inter_start=log(interpolationMax(zz)*x_max/scale - interpolationRange); 
        log_inter_end=log(interpolationMax(zz)*x_max/scale);

        for kk=1:numel(pos_out) 
            if log(pos_out(kk)) > log_inter_start;
                inter_start_index = kk;
                break;
            end      
        end

        for kk=1:numel(pos_out) 
            if log(pos_out(kk)) > log_inter_end;
                inter_end_index = kk;
                break;
            end      
        end
        
        inter_start_index = inter_end_index - 32;
        %fprintf('%d\n', inter_start_index - inter_end_index);
        
        xx = log(real(pos_out(inter_start_index:inter_end_index)));
        yy = real( log(real(correlationFn(inter_start_index:inter_end_index)) ) );

        [fitobject, gof] = fit(xx, yy, 'poly1');

        slope(zz) = fitobject.p1;
        intercept(zz) = fitobject.p2;
        sse(zz) = gof.sse;
        rmse(zz) = gof.rmse;
        rsquare(zz) = gof.rsquare;
        adjrsquare(zz) = gof.adjrsquare;
    
        confint95{zz} = confint(fitobject, 0.95);
        confint99{zz} = confint(fitobject, 0.9973);
    
        interpolationMidpoint(zz) = 0.5*( pos_out(inter_start_index) + pos_out(inter_end_index)  )/(x_max/scale);
    
    end
    fprintf('\n');
    
    confint99_upper = zeros(1, numel(interpolationMax));
    confint99_lower = zeros(1, numel(interpolationMax));
    
    for jj=1:numel(interpolationMax)
    
        confint99_lower(jj) = confint99{jj}(1,1);
        confint99_upper(jj) = confint99{jj}(2,1);

    end

    confint99_lower = slope - confint99_lower;
    confint99_upper = confint99_upper - slope;
    
    figureName = sprintf('exponent vs. scale');
    figure('Name', figureName);
    xlabel('$s$', 'Interpreter', 'LaTex'); 
    ylabel('$2 \Delta $', 'Interpreter', 'LaTex');

    plot(interpolationMax, slope, '-r');
    hold all;
    plot(interpolationMax, slope - confint99_lower, '-r');
    plot(interpolationMax, slope + confint99_upper, '-r');
    
    hold off;

end


