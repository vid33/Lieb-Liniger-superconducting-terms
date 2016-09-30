%CALCULATE EXPONENTS BY D-SCALING USING FULL EXPECTATION VALUE (NOT JUST
%DOMINANT CONTRIBUTIONS), USES MATLAB BUILT-IN ODE FUNCTIONS TO CALCULATE
%CORRELATION FUNCTIONS

clear;

uu = -5;
potential = 10;
interaction = 0;

Lambda = sqrt(20);
pi_prefactor = (-1i/2)*sqrt(2*Lambda);

plot_delta_vs_scale99 = true; %short for 99.73% conf interval
plot_delta_vs_scale95 = true;

field = false;
PI = false;
dPI = false; %first descendant of d_z Phi \equiv ddPhi
dddPhi = false; %second descendant of d_z Phi
dPhi = false;
hamiltonian = false;
density = false;
ekin = false;
ekin_density = false;
vertexOp = true;

%we use this eigenvalue of T to determine correlation lenght and do
%D-scaling
correlationEigenvalue = 2;

% for the vertex operator : exp(i beta \Phi) :
beta = (sqrt(4*pi));

%dim = 4:1:64;
dim = 4:1:64;


%exclude those that have not converged as much
%dim = [ 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26 , 27, 28, 29, 30, ...
%    31, 33, 34, 35, 36, 37, 38, 39, 41, 42, 43, 45, 46, 47, 49, 51, 53, 55, 57, 58, 59, 61, 62, 63 ];

eigval_start=29; %start interpolation at D=""+ 3

scale = 0.1:0.01:5;
scale = [0, scale]; %initial conditions for ode are defined at x=0;

position = zeros(numel(dim) - eigval_start, numel(scale));
correlationFn = zeros(numel(dim) - eigval_start, numel(scale));

%do interpolation starting for all ranges starting from eigval_start:dim(end) to eigval_start+scan_end:dim(end) 
scanDelta = 1;
%scanDelta = 3;

slope = zeros(scanDelta, numel(scale)-1);
intercept = zeros(scanDelta, numel(scale)-1);
sse = zeros(scanDelta, numel(scale)-1);
rmse = zeros(scanDelta, numel(scale)-1);
rsquare = zeros(scanDelta, numel(scale)-1);
adjrsquare = zeros(scanDelta, numel(scale)-1);
    
confint95 = cell(scanDelta, numel(scale)-1);
confint99 = cell(scanDelta, numel(scale)-1);

fIn= sprintf('data_LLuu/lr_eigenvectors/uu=%dv=%dc=%deigenvector=%d.mat',uu,  potential, interaction , correlationEigenvalue);
load(fIn, 'mu');

for kk=eigval_start:numel(dim)
    
    fprintf('%d ', dim(kk));
    if (rem(kk,10)==0)
        fprintf('\n');
    end 
    
    fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', dim(kk), dim(kk), uu, potential, interaction); 
    load(fIn, 'R', 'Q', 'D', 'matr', 'l', 'r', 'occupationDensity');
    Rkin = Q*R - R*Q;      
    
    if (hamiltonian == true)
        eDensity_rescaled = (1/Lambda)*trace(Rkin*matr*Rkin') ...
            - ((Lambda^2)/4)*( trace(R*R*matr) + trace(matr*R'*R') ) ...
            + ((Lambda^2)/2)*trace(R*matr*R') ;
    end
    
    pos_in = real(-1/mu(kk))*scale;
    position(kk-eigval_start+1, :) = pos_in;

    options=odeset('RelTol',1e-10, 'AbsTol', 1e-13);

    %[pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), [x_min x_max], reshape(R, D^2, 1), options );
    if (field == true)
        [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), pos_in, reshape(R, D^2, 1));
    elseif (PI ==true)
        [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), pos_in, ...
                                                    reshape(pi_prefactor*(R - R'), D^2, 1), options );
    elseif (dPI ==true)
        [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), pos_in, ...
                                                    reshape(pi_prefactor*(Rkin - Rkin'), D^2, 1) );
    elseif (dPhi ==true)
        [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), pos_in, ...
                                       reshape((1/sqrt(2*Lambda)*(Rkin +Rkin') + pi_prefactor*(R - R') ), D^2, 1), options);
    elseif (dddPhi == true)
        QR = Q*R - R*Q; QQR = Q*QR- QR*Q; RQR = R*QR- QR*R;
        QQQR = Q*QQR - QQR*Q;  RQQR = R*QQR - QQR*R; RRQR = R*RQR - RQR*R;
        
        [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), pos_in, ...
                                       reshape(pi_prefactor*(QQR - QQR' + R'*RQR - RQR'*R) ...
                                       +(1/sqrt(2*Lambda))*(QQQR + QQQR' +2*R'*RQQR + 2*RQQR'*R + R'*R'*RRQR + RRQR'*R*R ...
                                            + QR'*RQR + RQR'*QR)   , D^2, 1), options);
        
    elseif (hamiltonian == true)
        [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), pos_in, ...
                            reshape( (1/Lambda)*(Rkin')*Rkin ...
                            -((Lambda^2)/4)*(R*R + R'*R') + ((Lambda^2)/2)*R'*R, D^2, 1 ) );
    elseif (density  == true)
        [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), pos_in, reshape(R'*R, D^2, 1) );
    elseif (ekin == true)
        [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), pos_in, reshape(Rkin, D^2, 1), options );
    elseif (ekin_density == true)
        [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), pos_in, reshape(Rkin'*Rkin, D^2, 1), options );
    elseif (vertexOp == true)
        %compute (exp(beta*R\times 1...) acting on left and right
        %eigenvectors of T, and exp..T acting on the l*exp(beta*R\times1...) 
        
        options=odeset('RelTol',1e-10, 'AbsTol', 1e-16);
                    
        [pos_out1, out_corr1] = ode113(@(t,x) vertex_action(t, x, beta, Lambda, R, D), [0 1], ...
                reshape(eye(D), D^2, 1), options ) ;
                                
        matrix_left_vertex = transpose(reshape(transpose(out_corr1(end, :)), D, D));
        
        [pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), pos_in, ...
                                    reshape(transpose(matrix_left_vertex), D^2, 1), options);
        
        [pos_out_right, out_corr_right] = ode113(@(t,x) vertex_action_right(t, x, -beta, Lambda, R, D), [0 1], ...
                reshape(transpose(matr), D^2, 1), options ) ;
        matrix_right_vertex = transpose(reshape(transpose(out_corr_right(end, :)), D, D));
            
    end
        
    for mm=1:numel(scale)
        
        matrix_left = transpose(reshape(transpose(out_corr(mm, :)), D, D));
        
        if (field == true)
            correlationFn(kk-eigval_start+1, mm) = trace(matrix_left*conj(R)*transpose(matr));
        elseif (PI == true)
            correlationFn(kk-eigval_start+1, mm) = trace(matrix_left*pi_prefactor*( transpose(matr)*transpose(R) ...
                                            - conj(R)*transpose(matr) ));
        elseif (dPI == true)
            correlationFn(kk-eigval_start+1, mm) = trace(matrix_left*pi_prefactor*( transpose(matr)*transpose(Rkin) ...
                                            - conj(Rkin)*transpose(matr) ));
        elseif (dPhi == true)
            correlationFn(kk-eigval_start+1, mm) = trace(matrix_left*( 1/sqrt(2*Lambda)*( transpose(matr)*transpose(Rkin) ...
                                + conj(Rkin)*transpose(matr) ) ...
                                + pi_prefactor*( transpose(matr)*transpose(R) - conj(R)*transpose(matr) ) ));
        elseif (dddPhi == true)
            correlationFn(kk-eigval_start+1, mm) = trace(matrix_left*( pi_prefactor*(transpose(matr)*transpose(QQR) - conj(QQR)*transpose(matr) ...
                                                        + conj(R)*transpose(matr)*transpose(RQR) - conj(RQR)*transpose(matr)*transpose(R) ) ...
                                                        +1/sqrt(2*Lambda)*( transpose(matr)*transpose(QQQR) + conj(QQQR)*transpose(matr) ...
                                                        +2*conj(R)*transpose(matr)*transpose(RQQR) + 2*conj(RQQR)*transpose(matr)*transpose(R) ...
                                                        +conj(R*R)*transpose(matr)*transpose(RRQR) + conj(RRQR)*transpose(matr)*transpose(R*R) ...
                                                        +conj(QR)*transpose(matr)*transpose(RQR) + conj(RQR)*transpose(matr)*transpose(QR) ) ) );
        elseif (hamiltonian == true)
        correlationFn(kk-eigval_start+1, mm) = trace(matrix_left*( (1/Lambda)*conj(Rkin)*transpose(matr)*transpose(Rkin) ...
                               -((Lambda^2)/4)*( conj(R*R)*transpose(matr) + transpose(matr)*transpose(R*R) ) ...
                               +((Lambda^2)/2)*conj(R)*transpose(matr)*transpose(R) ) ) ...
                               - eDensity_rescaled^2;
        elseif (density == true)
            correlationFn(kk-eigval_start+1, mm) =  ...
                (1/(occupationDensity^2))*(trace(matrix_left*conj(R)*transpose(matr)*transpose(R)) - occupationDensity^2); 
        elseif (ekin == true)
            correlationFn(kk-eigval_start+1, mm) = ...
                trace(matrix_left*conj(Rkin)*transpose(matr)) - trace(Rkin*matr)- trace(matr*Rkin');
        elseif (ekin_density == true)
            correlationFn(kk-eigval_start+1, mm) = ...
                trace(matrix_left*conj(Rkin)*transpose(matr)*transpose(Rkin)) - trace(Rkin*matr*Rkin')^2;
        elseif (vertexOp == true)
            correlationFn(kk-eigval_start+1, mm) = ...
                trace(matrix_left*transpose(matrix_right_vertex));

        end
    end
end

fprintf('\n');

for jj=2:numel(scale)
    
    for zz=1:scanDelta
    
    xx = log(position( zz:end, jj));
    yy = log(abs( real(correlationFn( zz:end, jj)) ));

    [fitobject, gof] = fit(xx, yy, 'poly1');

    slope(zz, jj-1) = fitobject.p1;
    intercept(zz, jj-1) = fitobject.p2;
    sse(zz, jj-1) = gof.sse;
    rmse(zz, jj-1) = gof.rmse;
    rsquare(zz, jj-1) = gof.rsquare;
    adjrsquare(zz, jj-1) = gof.adjrsquare;
    
    confint95{zz, jj-1} = confint(fitobject, 0.95);
    confint99{zz, jj-1} = confint(fitobject, 0.9973);

end

end

confint99_upper = zeros(scanDelta, numel(scale)-1);
confint99_lower = zeros(scanDelta, numel(scale)-1);
for jj=2:numel(scale)

    for zz=1:scanDelta    
        confint99_lower(zz, jj-1) = confint99{zz, jj-1}(1,1);
        confint99_upper(zz, jj-1) = confint99{zz, jj-1}(2,1);
    end

end
confint99_lower = slope - confint99_lower;
confint99_upper = confint99_upper - slope;

confint95_upper = zeros(scanDelta, numel(scale)-1);
confint95_lower = zeros(scanDelta, numel(scale)-1);
for jj=2:numel(scale)
    for zz=1:scanDelta
        confint95_lower(zz, jj-1) = confint95{zz, jj-1}(1,1);
        confint95_upper(zz, jj-1) = confint95{zz, jj-1}(2,1);
    end
end
confint95_lower = slope - confint95_lower;
confint95_upper = confint95_upper - slope;

scale = scale(2:end);

if (plot_delta_vs_scale99 == true)

    figureName = sprintf('exponent vs. scale 99');
    figure('Name', figureName);
    hold all;
    xlabel('$s$', 'Interpreter', 'LaTex'); 
    ylabel('$2 \Delta $', 'Interpreter', 'LaTex');
    
    for zz=1:scanDelta

    plot(scale, slope(zz, :), '-r');
    plot(scale, slope(zz, :) - confint99_lower(zz, :), '-k');
    plot(scale, slope(zz, :) + confint99_upper(zz, :), '-k');

    
    end
    %errorbar( scale, slope, confint99_lower, confint99_upper, 'x', 'MarkerEdgeColor', 'k');
    
   % exactValue = -0.6366;
    exactValue = -1.432;
    %exactValue = -0.1668575;        
    line([scale(1) ; scale(end)],[exactValue ; exactValue], 'Color','g', 'LineWidth', 1);
    hold off;

end

if (plot_delta_vs_scale95 == true)

    figureName = sprintf('exponent vs. scale 95');
    figure('Name', figureName);
    hold all;
    xlabel('$s$', 'Interpreter', 'LaTex'); 
    ylabel('$2 \Delta $', 'Interpreter', 'LaTex');
    
    for zz=1:scanDelta

    plot(scale, slope(zz, :), '-r');
    plot(scale, slope(zz, :) - confint95_lower(zz, :), '-k');
    plot(scale, slope(zz, :) + confint95_upper(zz, :), '-k');

    
    end
    %errorbar( scale, slope, confint99_lower, confint99_upper, 'x', 'MarkerEdgeColor', 'k');
% exactValue = -0.6366;
    %exactValue = -4;
    exactValue = -1.432;
    %exactValue = -4;
    %exactValue = -2;
    %exactValue = -0.1668575;        
    line([scale(1) ; scale(end)],[exactValue ; exactValue], 'Color','g', 'LineWidth', 1);
    hold off;

end

