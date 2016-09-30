clear;

dimension = 32:64;

uu = -5;
potential = 10;
interaction = 0;

Lambda = sqrt(20);
pi_prefactor = (1i/2)*sqrt(2*Lambda);

%used to set maximum x value in terms of correlation length
scale = 4;

field = false;
PI = false; %not pi, to avoid conflict with 3.1415..
dPI = false;
dPhi = false;
hamiltonian = false;
density = false;
ekin = false;
ekin_density = false;
vertexOp = true;

% for the vertex operator : exp(i beta \Phi) :
if (vertexOp == true)
    beta = 1;
end
    
interpolationMax = 0.2:0.0075:2; %max here should be less than scale
interpolationMidpoint = zeros(numel(dimension), numel(interpolationMax));

slope = zeros(numel(dimension), numel(interpolationMax)); 
intercept = zeros(numel(dimension), numel(interpolationMax));
sse = zeros(numel(dimension), numel(interpolationMax)); 
rmse = zeros(numel(dimension), numel(interpolationMax));
rsquare = zeros(numel(dimension), numel(interpolationMax)); 
adjrsquare = zeros(numel(dimension), numel(interpolationMax));
confint95 = cell(numel(dimension), numel(interpolationMax)); 
confint99 = cell(numel(dimension), numel(interpolationMax));

confint99_upper = zeros(numel(dimension), numel(interpolationMax));
confint99_lower = zeros(numel(dimension), numel(interpolationMax));

x_max = zeros(1, numel(dimension));
x_min = 0;

%load correlation length values
for ff=1:numel(dimension)
    D = dimension(ff);
    fIn= sprintf('data_LLuu/lr_eigenvectors/uu=%dv=%dc=%deigenvector=%d.mat',uu,  potential, interaction , 3);
    load(fIn, 'mu', 'dim');

    for kk=1:numel(dim)
        if (dim(kk) == D)
            x_max(ff) = real(-1/mu(kk))*scale;
        end
    end
end

clearvars 'dim';

for ff=1:numel(dimension)
    
    D = dimension(ff);
    
    fprintf('D=%d\n', D);

    fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction); 
    load(fIn, 'R', 'Q', 'Rkin', 'matr');
    Rkin = Q*R - R*Q;

    %For hamiltonian - hamiltonian correlator
    if (hamiltonian == true)
        eDensity_rescaled = (1/Lambda)*trace(Rkin*matr*Rkin') ...
          - ((Lambda^2)/4)*( trace(R*R*matr) + trace(matr*R'*R') ) ...
          + ((Lambda^2)/2)*trace(R*matr*R') ;

    end

    options=odeset('RelTol',1e-10, 'AbsTol', 1e-13);

    if (field ==true)
        [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), [x_min x_max(ff)], reshape(R, D^2, 1), options );
    elseif (PI ==true)
        [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), [x_min x_max(ff)], ...
                                                    reshape(pi_prefactor*(R - R'), D^2, 1), options );
    elseif (dPI ==true)
        [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), [x_min x_max(ff)], ...
                                                    reshape(pi_prefactor*(Rkin - Rkin'), D^2, 1), options );
    elseif (dPhi ==true)
        [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), [x_min x_max(ff)], ...
                                       reshape( 0.5*(1/sqrt(2*Lambda)*(Rkin + Rkin') + pi_prefactor*(R - R') ), D^2, 1), options );
    elseif (hamiltonian == true)
        [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), [x_min x_max(ff)], ...
                reshape( (1/Lambda)*(Rkin')*Rkin ...
                -((Lambda^2)/4)*(R*R + R'*R') + ((Lambda^2)/2)*R'*R, D^2, 1 ) );
    elseif (density == true)
        [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), [x_min x_max(ff)], reshape(R'*R, D^2, 1), options );
    elseif (ekin == true)
        [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), [x_min x_max(ff)], reshape(Rkin, D^2, 1), options );
    elseif (ekin_density == true)
        [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), [x_min x_max(ff)], reshape(Rkin'*Rkin, D^2, 1), options );
    elseif (vertexOp == true)
        %compute (exp(beta*R\times 1...) acting on left and right
        %eigenvectors of T, and exp..T acting on the l*exp(beta*R\times1...) 
                    
        [pos_out1, out_corr1] = ode113(@(t,x) vertex_action(t, x, beta, Lambda, R, D), [0 1], ...
                reshape(eye(D), D^2, 1), options ) ;
                                
        matrix_left_vertex = transpose(reshape(transpose(out_corr1(end, :)), D, D));
        
        [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), [x_min x_max(ff)], ...
                                    reshape(transpose(matrix_left_vertex), D^2, 1), options);
        
        [pos_out_right, out_corr_right] = ode113(@(t,x) vertex_action_right(t, x, -beta, Lambda, R, D), [0 1], ...
                reshape(transpose(matr), D^2, 1), options ) ;
        matrix_right_vertex = transpose(reshape(transpose(out_corr_right(end, :)), D, D));
    end
    
    correlationFn = zeros(1, numel(pos_out));

    for kk=1:numel(pos_out)
        matrix_left = transpose(reshape(transpose(out_corr(kk, :)), D, D));
    
        if (field == true)
            correlationFn(kk) = trace(matrix_left*conj(R)*transpose(matr)) - trace(R*matr)*trace(matr*R');
        elseif (PI == true)
            correlationFn(kk) = trace(matrix_left*pi_prefactor*( transpose(matr)*transpose(R) - conj(R)*transpose(matr) ));
        elseif (dPI == true)
            correlationFn(kk) = trace(matrix_left*pi_prefactor*( transpose(matr)*transpose(Rkin) ...
                        - conj(Rkin)*transpose(matr) ));
        elseif (dPhi == true)
            correlationFn(kk) = trace(matrix_left*0.5*( 1/sqrt(2*Lambda)*( transpose(matr)*transpose(Rkin) + conj(Rkin)*transpose(matr) ) ...
                + pi_prefactor*( transpose(matr)*transpose(R) - conj(R)*transpose(matr) ) ) );
        elseif (hamiltonian == true)
            correlationFn(kk) = trace(matrix_left*( (1/Lambda)*conj(Rkin)*transpose(matr)*transpose(Rkin) ...
                               -((Lambda^2)/4)*( conj(R*R)*transpose(matr) + transpose(matr)*transpose(R*R) ) ...
                               +((Lambda^2)/2)*conj(R)*transpose(matr)*transpose(R) ) ) ...
                               - eDensity_rescaled^2;
        elseif (density == true)
            correlationFn(kk) = ...
                (1/(occupationDensity^2))*(trace(matrix_left*conj(R)*transpose(matr)*transpose(R)) - occupationDensity^2);
        elseif (ekin == true)
            correlationFn(kk) = trace(matrix_left*conj(Rkin)*transpose(matr)) - trace(Rkin*matr)- trace(matr*Rkin');
        elseif (ekin_density == true)
         correlationFn(kk) = ...
            trace(matrix_left*conj(Rkin)*transpose(matr)*transpose(Rkin)) - trace(Rkin*matr*Rkin')^2;
        elseif (vertexOp == true)
            correlationFn(kk) = ...
                trace(matrix_left*transpose(matrix_right_vertex));
        end
        
    end


    correlationFn = transpose(correlationFn);
    
    rangeScale = 0.01;
    interpolationRange = rangeScale*x_max(ff)/scale;
    
    fprintf('Performing interpolations \n');
    for zz=1:numel(interpolationMax)
        fprintf('%d ', zz);
        if (rem(zz,10)==0)
           fprintf('\n');
        end
        
        log_inter_start=log(interpolationMax(zz)*x_max(ff)/scale - interpolationRange); 
        log_inter_end=log(interpolationMax(zz)*x_max(ff)/scale);

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

        slope(ff, zz) = fitobject.p1;
        intercept(ff, zz) = fitobject.p2;
        sse(ff, zz) = gof.sse;
        rmse(ff, zz) = gof.rmse;
        rsquare(ff, zz) = gof.rsquare;
        adjrsquare(ff, zz) = gof.adjrsquare;
    
        confint95{ff, zz} = confint(fitobject, 0.95);
        confint99{ff, zz} = confint(fitobject, 0.9973);
    
        interpolationMidpoint(ff, zz) = 0.5*( pos_out(inter_start_index) + pos_out(inter_end_index)  )/(x_max(ff)/scale);
    
    end
    fprintf('\n');
    
    for jj=1:numel(interpolationMax) 
        confint99_lower(ff, jj) = confint99{ff, jj}(1,1);
        confint99_upper(ff, jj) = confint99{ff, jj}(2,1);
    end
    
end


figureName = sprintf('exponent vs. scale');
figure('Name', figureName);
xlabel('$s$', 'Interpreter', 'LaTex'); 
ylabel('$2 \Delta $', 'Interpreter', 'LaTex');

hold all;
    
confint99_lower = slope - confint99_lower;
confint99_upper = confint99_upper - slope;
    
for yy = 1: numel(dimension)
    
    plot(interpolationMax, slope(yy, :), '-k');
    
    %plot(interpolationMax, slope - confint99_lower, '-r');
    %plot(interpolationMax, slope + confint99_upper, '-r');
    
end
    
hold off;


figureName = sprintf('exponent vs. position');
figure('Name', figureName);
xlabel('pos', 'Interpreter', 'LaTex'); 
ylabel('$2 \Delta $', 'Interpreter', 'LaTex');

hold all;
    
confint99_lower = slope - confint99_lower;
confint99_upper = confint99_upper - slope;
    
for yy = 1: numel(dimension)
    
    plot(interpolationMax*x_max(yy)/scale, slope(yy, :), '-k');
    
    %plot(interpolationMax, slope - confint99_lower, '-r');
    %plot(interpolationMax, slope + confint99_upper, '-r');
    
end
    
hold off;


