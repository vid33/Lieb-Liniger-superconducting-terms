clear;

uu = -5;
potential = 10;
interaction = 0;
%eigenvectors = [45 ]; 

Lambda = sqrt(20);
pi_prefactor = (1i/2)*sqrt(2*Lambda);

PI = false;
dPI = false; % == ddPhi
dPhi = false;
dddPhi = false;
hamiltonian = false;
vertex = false;
dPhi_vertex = false;
dPhi_hamiltonian = true;

% for the vertex operator : exp(i beta \Phi) :
beta = 1;

eigenvectorNo_max = 500;
fIn = sprintf('data_LLuu/lr_eigenvectors/%dEIGVECSuu=%dv=%dc=%d.mat',eigenvectorNo_max, uu,  potential, interaction);
load(fIn,'mu', 'matl2', 'matr2');

%mu = cell(eigenvectorNo_max, numel(dim));

DIM = 16;

mu_2 = mu{2, DIM-1};

eigenvectorNo_max = min(500, DIM^2);

contribution = zeros(1, eigenvectorNo_max);
contribution_l = zeros(1, eigenvectorNo_max);
contribution_r = zeros(1, eigenvectorNo_max);

fIn_tmp = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', DIM, DIM, uu, potential, interaction); 
load(fIn_tmp, 'R', 'Q', 'D', 'matr');
Rkin = Q*R - R*Q;

contributionTotal = 0;
            
for pp=1:eigenvectorNo_max
    
    fprintf('%d ', pp);
    if (rem(pp,10)==0)
                fprintf('\n');
    end 
    
    if PI == true
        tmp1 =  trace(R*matr2{pp, DIM-1}) -trace(matr2{pp, DIM-1}*R') ; 
        tmp2 = trace(matr*transpose(matl2{pp, DIM-1})*R) ...
                - trace(matr*R'*transpose(matl2{pp, DIM-1}));

    elseif dPhi == true
        tmp1 = 0.5*( 1/sqrt(2*Lambda)*( trace(Rkin*matr2{pp, DIM-1}) + trace(matr2{pp, DIM-1}*Rkin' ) )  ...
                    + pi_prefactor*( trace(R*matr2{pp, DIM-1}) - trace(matr2{pp, DIM-1}*R') ) );
        tmp2 = 0.5*( 1/sqrt(2*Lambda)*( trace(matr*transpose(matl2{pp, DIM-1})*Rkin) + trace(matr*Rkin'*transpose(matl2{pp, DIM-1})) ) ...
                    + pi_prefactor*( trace(matr*transpose(matl2{pp, DIM-1})*R) - trace(matr*R'*transpose(matl2{pp, DIM-1})) ) );
     
    elseif dPI == true
         %[pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), pos_in, ...
          %                                          reshape(pi_prefactor*(Rkin - Rkin'), D^2, 1) );
        tmp1 =  trace(Rkin*matr2{pp, DIM-1}) -trace(matr2{pp, DIM-1}*Rkin') ; 
        tmp2 = trace(matr*transpose(matl2{pp, DIM-1})*Rkin) ...
                - trace(matr*Rkin'*transpose(matl2{pp, DIM-1}));
    elseif dddPhi == true
        
        %NOT COMPLETED
         QR = Q*R - R*Q; QQR = Q*QR- QR*Q; RQR = R*QR- QR*R;
        QQQR = Q*QQR - QQR*Q;  RQQR = R*QQR - QQR*R; RRQR = R*RQR - RQR*R;
        
        %[pos_out, out_corr] = ode113(@(t,x) T_action_2(t, x, R, Q, D), pos_in, ...
         %                              reshape(pi_prefactor*(QQR - QQR' + R'*RQR - RQR'*R) ...
         %                              +(1/sqrt(2*Lambda))*(QQQR + QQQR' +2*R'*RQQR + 2*RQQR'*R + R'*R'*RRQR + RRQR'*R*R ...
         %                                   + QR'*RQR + RQR'*QR)   , D^2, 1), options);
         
        %tmp1 = 
         
         
    elseif hamiltonian == true
        tmp1 = (1/Lambda)*( trace(Rkin*matr2{pp, DIM-1}*Rkin') - ((Lambda^2)/4)*( trace(R*R*matr2{pp, DIM-1}) ...
            + trace(matr2{pp, DIM-1}*R'*R') ) + ((Lambda^2)/2)*trace(R*matr2{pp, DIM-1}*R')  );
                    %+ 1i*(1/2)*sqrt(2*Lambda)*( ...
                    %    - trace(Rkin*matr2{k}*R') + trace(R*matr2{k}*Rkin') ) );
                    % TOT DERIVATIVE             + trace(Rkin2*matr2{k})  - trace(matr2{k}*Rkin2')) );
        tmp2=      (1/Lambda)*( trace(matr*Rkin'*transpose(matl2{pp, DIM-1})*Rkin)   ...
                    - ((Lambda^2)/4)*( trace(matr*transpose(matl2{pp, DIM-1})*R*R)  + trace(matr*R'*R'*transpose(matl2{pp, DIM-1})) ) ...
                    + ((Lambda^2)/2)*trace(matr*R'*transpose(matl2{pp, DIM-1})*R) );
                    %+ 1i*(1/2)*sqrt(2*Lambda)*( ...
                    %          - trace(matr*R'*transpose(matl2{k})*Rkin) + trace(matr*Rkin'*transpose(matl2{k})*R) ) );
                    % TOT. DERIVATIVES     + trace(matr*transpose(matl2{k})*Rkin2) - trace(matr*Rkin'*transpose(matl2{k}))    ) );
       % contribution(pp) = exp(-1*mu{pp, DIM-1}/mu_2 )*tmp1*tmp2;
    elseif vertex == true
        options=odeset('RelTol',1e-10, 'AbsTol', 1e-13);
                    
        [pos_out, out_corr] = ode113(@(t,x) vertex_action(t, x, beta, Lambda, R, D), [0 1], ...
                reshape(eye(D), D^2, 1), options ) ;
                                
        matrix_left = transpose(reshape(transpose(out_corr(end, :)), D, D));
        
        [pos_out_right, out_corr_right] = ode113(@(t,x) vertex_action_right(t, x, -beta, Lambda, R, D), [0 1], ...
                reshape(transpose(matr), D^2, 1), options ) ;
        matrix_right_vertex = transpose(reshape(transpose(out_corr_right(end, :)), D, D));
                
        tmp1 = trace(transpose(matrix_left)*matr2{pp, DIM-1});
        tmp2 = trace(transpose(matl2{pp, DIM-1})*matrix_right_vertex);
    elseif dPhi_vertex == true
        options=odeset('RelTol',1e-10, 'AbsTol', 1e-13);
        tmp1 = 0.5*( 1/sqrt(2*Lambda)*( trace(Rkin*matr2{pp, DIM-1}) + trace(matr2{pp, DIM-1}*Rkin' ) )  ...
            + pi_prefactor*( trace(R*matr2{pp, DIM-1}) - trace(matr2{pp, DIM-1}*R') ) );
        
        [pos_out_right, out_corr_right] = ode113(@(t,x) vertex_action_right(t, x, -beta, Lambda, R, D), [0 1], ...
                reshape(transpose(matr), D^2, 1), options ) ;
        matrix_right_vertex = transpose(reshape(transpose(out_corr_right(end, :)), D, D));
        tmp2 = trace(transpose(matl2{pp, DIM-1})*matrix_right_vertex);
    elseif dPhi_hamiltonian == true
        tmp1 = 0.5*( 1/sqrt(2*Lambda)*( trace(Rkin*matr2{pp, DIM-1}) + trace(matr2{pp, DIM-1}*Rkin' ) )  ...
            + pi_prefactor*( trace(R*matr2{pp, DIM-1}) - trace(matr2{pp, DIM-1}*R') ) );
        tmp2=      (1/Lambda)*( trace(matr*Rkin'*transpose(matl2{pp, DIM-1})*Rkin)   ...
            - ((Lambda^2)/4)*( trace(matr*transpose(matl2{pp, DIM-1})*R*R)  + trace(matr*R'*R'*transpose(matl2{pp, DIM-1})) ) ...
            + ((Lambda^2)/2)*trace(matr*R'*transpose(matl2{pp, DIM-1})*R) );

    end     
    
    contribution_l(pp) = tmp1;
    contribution_r(pp) = tmp2;
    contribution(pp) = tmp1*tmp2;
  %  contribution(pp) = exp(-1*mu{pp, DIM-1}/mu_2 )*tmp1*tmp2;
                
end

plot_x = 1:1:eigenvectorNo_max;

figure; plot(plot_x, (contribution), 'x');
figure; plot(plot_x, contribution_l, 'xb', plot_x, contribution_r, 'xr');
