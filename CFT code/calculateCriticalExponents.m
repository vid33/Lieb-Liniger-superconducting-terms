clear;

uu = -5;
potential = 10;
interaction = 0;

eigenvectors = [4 ]; %nth eigenvector
correlationEigenvalue = 2;  %use this eigenvalue to determine correlation length
scale=0.1:0.05:5;  %range of scales at which we calculate the correlator
%scale = 1;

eigval_start = 29; % position in dim where we start the interpolation

%if we wanna end befor maximum available dimension = dim(end)
offset = 0;

generatePlot = true;
plot_delta_vs_scale = true;

%%%%%%%%%%%%%%%%%%%%%%%%%% SELECT ONE!
field = false;
density = false;
ekin = false;
ekinDensity = false;
uuTerm = false;
jacobi = false;

phi = false;
PI = true;
DPI = false;
DPhi = false;
DDPhi = false;

TDPhi = false;
TDDPhi = false;

hamiltonian = false;

vertexOp = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%

kappa = 6/(sqrt(12)+1); 
Lambda = sqrt(20);

% for the vertex operator : exp(i beta \Phi) :
beta = 3.5;

%fIn = sprintf('data_LLuu/lr_eigenvectors/uu=%dv=%dc=%deigenvector=%d.mat',uu,  potential, interaction , 2);
eigenvectorNo_max = 15;  
fIn = sprintf('data_LLuu/lr_eigenvectors/%dEIGVECSuu=%dv=%dc=%d.mat',eigenvectorNo_max, uu,  potential, interaction);
load(fIn, 'dim', 'mu', 'matl2', 'matr2');

dim = dim(eigval_start:end-offset);

%Set eigenvalue used for correlation length
mu_2 = zeros(1, numel(dim));
for kk=1:numel(dim)
    mu_2(kk) = mu{correlationEigenvalue,eigval_start+kk-1};
end

dim_lowest = dim(1); fprintf('DIM lowest %d\n', dim_lowest);


%Disconnected contributions

aal_disconn = zeros(1, numel(dim));
aar_disconn = zeros(1, numel(dim));


for k=1:numel(dim)
    
    fIn_tmp = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', dim(k), dim(k), uu, potential, interaction); 
    load(fIn_tmp, 'R', 'Q', 'D', 'matr', 'l', 'r', 'occupationDensity');
    Rkin = Q*R - R*Q;
    
    if (field == true)
        aal_disconn(k) = trace(matr*R');
        aar_disconn(k) = trace(R*matr);
   
    elseif (density == true)
        aal_disconn(k) = trace(R*matr*R');
        aar_disconn(k) = trace(matr*R'*R);
        
    elseif (ekin == true)   
        aal_disconn(k) = trace(matr*Rkin');
        aar_disconn(k) = trace(matr*Rkin);
        
    elseif (ekinDensity == true)
        Rkin = Q*R - R*Q;
        
        aal_disconn(k) = trace(Rkin*matr*Rkin');
        aar_disconn(k) = trace(matr*(Rkin')*Rkin);
   
    elseif (phi == true)
        aal_disconn(k) = trace(R*matr) + trace(matr*R');
        aar_disconn(k) = trace(matr*R) + trace(matr*R');
   
    elseif (PI == true)
        aal_disconn(k) = trace(R*matr) -trace(matr*R');
        aar_disconn(k) = trace(matr*R) - trace(matr*R');
        
    elseif (DPI == true)
        aal_disconn(k) = trace(Rkin*matr) -trace(matr*Rkin');
        aar_disconn(k) = trace(matr*Rkin) - trace(matr*Rkin');
    
    elseif (DPhi == true)
        DDR = Q*Rkin - Rkin*Q;
    
        aal_disconn(k) = 0.5*( 1/sqrt(2*sqrt(Lambda))*( trace(Rkin*matr) + trace(matr*Rkin') )...
            + 0.5*1i*sqrt(2*sqrt(Lambda))*( trace(R*matr) - trace(matr*R') ) ) ;
        aar_disconn(k) = 0.5*(  1/sqrt(2*sqrt(Lambda))*( trace(matr*Rkin) - trace(matr*Rkin') ) ...
            + 0.5*1i*sqrt(2*sqrt(Lambda))*( trace(matr*R) - trace(matr*R') )  );    
    
    elseif (hamiltonian == true)
        aal_disconn(k) = (1/Lambda)*( trace(Rkin*matr*Rkin') ...
            - ((Lambda^2)/4)*( trace(R*R*matr) + trace(matr*R'*R') ) ...
            + ((Lambda^2)/2)*trace(R*matr*R') ); 
            %    + 4i*Lambda*( - trace(Rkin*matr*R') + trace(R*matr*Rkin') );
        aar_disconn(k) = (1/Lambda)*(  trace(matr*(Rkin')*Rkin)   ...
            - ((Lambda^2)/4)*( trace(matr*R*R)  + trace(matr*R'*R') ) ...
            + ((Lambda^2)/2)*trace(matr*R'*R)  );
            %+ 4i*Lambda*( ...
            %        - trace(matr*R'*Rkin) + trace(matr*Rkin'*R) );
    
    elseif (DDPhi == true)      
        aal_disconn(k) = 0.5*1i*sqrt(2*Lambda)*( trace(Rkin*matr) - trace(matr*Rkin') );
        aar_disconn(k) = 0.5*1i*sqrt(2*Lambda)*( trace(matr*Rkin) - trace(matr*Rkin') );
    
    
    elseif(vertexOp == true)
        fprintf('Current dim %d\n', dim(k));
                    
        %aal_disconn(k) =l*expm( (-beta/sqrt(2*sqrt(20)) )*1i*( kron(R, eye(D)) + kron(eye(D), conj(R))))*r;             
        %aar_disconn(k) = l*expm(( beta/sqrt(2*sqrt(20)) )*1i*( kron(R, eye(D)) + kron(eye(D), conj(R))))*r;
        
        options=odeset('RelTol',1e-10, 'AbsTol', 1e-13);
                    
        [pos_out, out_corr] = ode113(@(t,x) vertex_action(t, x, beta, Lambda, R, D), [0 1], ...
                reshape(eye(D), D^2, 1), options ) ;
                                
        matrix_left = transpose(reshape(transpose(out_corr(end, :)), D, D));
                
        aal_disconn(k) = trace(transpose(matrix_left)*matr);
            
        [pos_out, out_corr] = ode113(@(t,x) vertex_action(t, x, -beta, Lambda, R, D), [0 1], ...
                reshape(eye(D), D^2, 1), options ) ;
            
        matrix_left = transpose(reshape(transpose(out_corr(end, :)), D, D));
            
        aar_disconn(k) = trace(transpose(matrix_left)*matr);
        
    end
end

   
%%Calculate overlaps with pp-th left and right eigenvectors

aal = zeros(numel(eigenvectors), numel(dim));
aar = zeros(numel(eigenvectors), numel(dim));

for pp=1:numel(eigenvectors)
    
    fprintf('EIGENVECTOR %d\n', eigenvectors(pp));
    for k=1:numel(dim)
        
        fprintf('%d ', k);
        if (rem(k,10)==0)
                fprintf('\n');
        end 
            
        dimIndex = k-1+eigval_start;
        eigv = eigenvectors(pp);
    
        fIn_tmp = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', dim(k), dim(k), uu, potential, interaction); 
        load(fIn_tmp, 'R', 'Q', 'D', 'matr', 'occupationDensity');   
        Rkin = Q*R - R*Q;
                   
        if (field == true)
            aal(pp, k) = aal(pp, k) + exp(bb*0.5*real(-1*mu{pp, k}/mu_2(k)))*trace(matr2{pp, k}*R');
            aar(pp, k) = aar(pp, k) + exp(bb*0.5*real(-1*mu{pp, k}/mu_2(k)))*trace(matr*transpose(matl2{pp, k})*R);
            
        elseif (density == true)
            tmp1 = trace(R*matr2{eigv, dimIndex}*R');
            tmp2 = trace(matr*R'*transpose(matl2{eigv, dimIndex})*R);
    
        elseif (ekin == true)    
            tmp1 = 0.5*( 1/sqrt(2*Lambda)*( trace(Rkin*matr2{eigv, dimIndex}) + trace(matr2{eigv, dimIndex}*Rkin' ) ) );
            tmp2 = 0.5*( 1/sqrt(2*Lambda)*( trace(matr*transpose(matl2{eigv, dimIndex})*Rkin)...
                + trace(matr*Rkin'*transpose(matl2{eigv, dimIndex})) ) );
                
        elseif (ekinDensity == true)   
            tmp1 = trace(Rkin*matr2{eigv, dimIndex}*Rkin');
            tmp2 = trace(matr*Rkin'*transpose(matl2{eigv, dimIndex})*Rkin);
                
        elseif (uuTerm == true)
            tmp1 = trace(R*R*matr2{eigv, dimIndex}) + trace(matr2{eigv, dimIndex}*R'*R');
            tmp2 = trace(matr*transpose(matl2{eigv, dimIndex})*R*R)  + trace(matr*R'*R'*transpose(matl2{eigv, dimIndex}));
        
        elseif (phi == true)
            aal(pp,k) = trace(R*matr2{k}) + trace(matr2{k}*R') ;
            aar(pp,k) = trace(matr*transpose(matl2{k})*R) + trace(matr*R'*transpose(matl2{k}));
        
        elseif (PI == true)
            aal(pp,k) =  trace(R*matr2{eigv, dimIndex}) - trace(matr2{eigv, dimIndex}*R') ; 
            aar(pp,k) = trace(matr*transpose(matl2{eigv, dimIndex})*R) - trace(matr*R'*transpose(matl2{eigv, dimIndex}));
            
        elseif (DPI == true)
            aal(pp,k) =  trace(Rkin*matr2{eigv, dimIndex}) - trace(matr2{eigv, dimIndex}*Rkin') ; 
            aar(pp,k) = trace(matr*transpose(matl2{eigv, dimIndex})*Rkin) - trace(matr*Rkin'*transpose(matl2{eigv, dimIndex}));
                   
        elseif (DPhi == true)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% D phi = d/dx ( Psi + Psi^dag) +
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% i*sqrt(20)(Psi = Psi^dag)

            aal(pp,k) = 0.5*( 1/sqrt(2*Lambda)*( trace(Rkin*matr2{eigv, dimIndex}) + trace(matr2{eigv, dimIndex}*Rkin' ) ) ...
                + 0.5*1i*sqrt(2*Lambda)*( trace(R*matr2{eigv, dimIndex}) - trace(matr2{eigv, dimIndex}*R') ) );
            aar(pp,k) = 0.5*( 1/sqrt(2*Lambda)*( trace(matr*transpose(matl2{eigv, dimIndex})*Rkin) + trace(matr*Rkin'*transpose(matl2{eigv, dimIndex})) ) ...
                + 0.5*1i*sqrt(2*Lambda)*( trace(matr*transpose(matl2{eigv, dimIndex})*R) - trace(matr*R'*transpose(matl2{eigv, dimIndex})) ) );
            
        elseif (DDPhi == true)  
            %DDR = Q*Rkin - Rkin*Q;
      
            aal(pp,k) = 0.5*1i*sqrt(2*Lambda)*( trace(Rkin*matr2{k}) - trace(matr2{k}*Rkin') );
            aar(pp,k) = 0.5*1i*sqrt(2*Lambda)*( trace(matr*transpose(matl2{k})*Rkin) - trace(matr*Rkin'*transpose(matl2{k})) );
      
        elseif (hamiltonian == true)
            aal(pp,k) = (1/Lambda)*( trace(Rkin*matr2{eigv, dimIndex}*Rkin') - ((Lambda^2)/4)*( trace(R*R*matr2{eigv, dimIndex}) ...
                    + trace(matr2{eigv, dimIndex}*R'*R') ) + ((Lambda^2)/2)*trace(R*matr2{eigv, dimIndex}*R')  );
                    %+ 1i*(1/2)*sqrt(2*Lambda)*( ...
                    %    - trace(Rkin*matr2{k}*R') + trace(R*matr2{k}*Rkin') ) );
                    % TOT DERIVATIVE             + trace(Rkin2*matr2{k})  - trace(matr2{k}*Rkin2')) );
            aar(pp,k)=      (1/Lambda)*( trace(matr*Rkin'*transpose(matl2{eigv, dimIndex})*Rkin)   ...
                    - ((Lambda^2)/4)*( trace(matr*transpose(matl2{eigv, dimIndex})*R*R)  + trace(matr*R'*R'*transpose(matl2{eigv, dimIndex})) ) ...
                    + ((Lambda^2)/2)*trace(matr*R'*transpose(matl2{eigv, dimIndex})*R) );
                    %+ 1i*(1/2)*sqrt(2*Lambda)*( ...
                    %          - trace(matr*R'*transpose(matl2{k})*Rkin) + trace(matr*Rkin'*transpose(matl2{k})*R) ) );
                    % TOT. DERIVATIVES     + trace(matr*transpose(matl2{k})*Rkin2) - trace(matr*Rkin'*transpose(matl2{k}))    ) );
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif(vertexOp == true)
            %aal_exp(k) = aal_exp(k) + exp(bb*0.5*real(-1*mu(k)/mu_2(k)))*(...
            %            l*expm((-sqrt(2*pi))*1i*beta*( kron(R, eye(D)) + kron(eye(D), conj(R))))*r2{k});
                    
            %aar_exp(k) = aar_exp(k) + exp(bb*0.5*real(-1*mu(k)/mu_2(k)))*(...
            %            l2{k}*expm(sqrt(2*pi)*1i*beta*( kron(R, eye(D)) + kron(eye(D), conj(R))))*r);
        
           % aal(pp, k) = l*expm( (-beta/sqrt(2*Lambda)  )*1i*( kron(R, eye(D)) + kron(eye(D), conj(R))))*r2{k};
                    
          %  aar(pp, k) = l2{k}*expm( ( beta/sqrt(2*Lambda) )*1i*( kron(R, eye(D)) + kron(eye(D), conj(R))))*r;
                   
            options=odeset('RelTol',1e-10, 'AbsTol', 1e-16);
                    
            [pos_out, out_corr] = ode113(@(t,x) vertex_action(t, x, beta, Lambda, R, D), [0 1], ...
                reshape(eye(D), D^2, 1), options ) ;
                                
            matrix_left = transpose(reshape(transpose(out_corr(end, :)), D, D));
                
            aal(pp, k) = trace(transpose(matrix_left)*matr2{eigv, dimIndex});
            
            [pos_out, out_corr] = ode113(@(t,x) vertex_action(t, x, -beta, Lambda, R, D), [0 1], ...
                reshape(transpose(matl2{eigv, dimIndex}), D^2, 1), options ) ;
            
            matrix_left = transpose(reshape(transpose(out_corr(end, :)), D, D));
            
            aar(pp, k) = trace(transpose(matrix_left)*matr);
                    
        elseif (TDPhi == true)
            aal(pp,k) = (Lambda/2)*( (1/Lambda)*trace(Rkin*matr2{k}*Rkin') ...
                        - (Lambda/4)*( trace(R*R*matr2{k}) + trace(matr2{k}*R'*R') ) + (Lambda/2)*trace(R*matr2{k}*R')  ...
                    + 1i*(1/2)*sqrt(2*Lambda)*( - trace(Rkin*matr2{k}*R') + trace(R*matr2{k}*Rkin') ) );
            aar(pp,k) = (  trace(matr*transpose(matl2{k})*Rkin) - trace(matr*Rkin'*transpose(matl2{k})) ...
                    + 1i*sqrt(20)*( trace(matr*transpose(matl2{k})*R) - trace(matr*R'*transpose(matl2{k})) ) ) ;
        
        elseif (jacobi == true)
            Rjacobi = R*Rkin - Rkin*R;
            aal(pp,k) = trace(matr2{k}*Rjacobi');
            aar(pp,k) = trace(matr*transpose(matl2{k})*Rjacobi);
    
                %aal_DDR(k) = trace(matr2{k}*DDR');
                %aar_DDR(k) = trace(matr*transpose(matl2{k})*DDR);
        end
        
    end
    
    fprintf('\n');
    
end

fprintf('\n');

%Calculate contributions at different correlation lengths

fitMeasure = zeros(1, numel(scale));
slope = zeros(1, numel(scale));
intercept = zeros(1, numel(scale));
sse = zeros(1, numel(scale));
rmse = zeros(1, numel(scale));
rsquare = zeros(1, numel(scale));
adjrsquare = zeros(1, numel(scale));
confint95 = cell(1, numel(scale));
confint99 = cell(1, numel(scale));
   
for zz=1:numel(scale)
     
    aa = zeros(1, numel(dim));
    
    fprintf('%d ', zz);
    if (rem(zz,10)==0)
                fprintf('\n');
    end 
     
    for pp=1:numel(eigenvectors)
        for k=1:numel(dim)
            
            dimIndex = k-1+eigval_start;
            eigv = eigenvectors(pp);
     
            aa(k) = aa(k) + exp(scale(zz)*real(-1*mu{eigv, dimIndex}/mu_2(k)))*aal(pp,k)*aar(pp,k);
   
        end
    end
    
    %xx = real(log(-1./mu));
    xx = real(log(-scale(zz)./mu_2));
    %yy = real(log(aal_disconn.*aar_disconn + aa));
    yy = real(log(aa));

    [fitobject, gof] = fit(transpose(xx), transpose(yy), 'poly1');

    slope(zz) = fitobject.p1;
    intercept(zz) = fitobject.p2;
    sse(zz) = gof.sse;
    rmse(zz) = gof.rmse;
    rsquare(zz) = gof.rsquare;
    adjrsquare(zz) = gof.adjrsquare;
    
    confint95{zz} = confint(fitobject, 0.95);
    confint99{zz} = confint(fitobject, 0.9973);
    
end
fprintf('\n');

[bestSSE, bestSSEIndex] = min(sse);
[bestRsquare, bestRsquareIndex] = max(rsquare);

bestScale = scale(bestSSEIndex);
bestSlope = slope(bestSSEIndex);
bestIntercept = intercept(bestSSEIndex);
bestConfInt95 = confint95{bestSSEIndex};
bestConfInt99 = confint99{bestSSEIndex};

fprintf('\n');


%%%PLOT STUFF

if (generatePlot == true)
    
    figureName = sprintf('a <-> field, disconnected term');
    figure('Name', figureName); plot(log(-1./mu_2),...
    log(aal_disconn.*aar_disconn), 'x');
    xlabel('$\log(\mu)$', 'Interpreter', 'LaTex'); 
    ylabel('$\log(a)$', 'Interpreter', 'LaTex');

    figureName = sprintf('a <-> field');
    figure('Name', figureName); plot(log(-1./mu_2), ...
            log(aa), 'x');
    xlabel('$\log(\mu)$', 'Interpreter', 'LaTex'); 
    ylabel('$\log(a)$', 'Interpreter', 'LaTex');
    
    figureName = sprintf('a <-> field, disconn + conn');
    figure('Name', figureName); plot(log(-1./mu_2), ...
        log(aa + aal_disconn.*aar_disconn), 'x');
    xlabel('$\log(\mu)$', 'Interpreter', 'LaTex'); 
    ylabel('$\log(a)$', 'Interpreter', 'LaTex');

end

confint99_upper = zeros(1, numel(scale));
confint99_lower = zeros(1, numel(scale));

for jj=1:numel(scale)
    
    confint99_lower(jj) = confint99{jj}(1,1);
    confint99_upper(jj) = confint99{jj}(2,1);

end

confint99_lower = slope - confint99_lower;
confint99_upper = confint99_upper - slope;

if (plot_delta_vs_scale == true)

    figureName = sprintf('exponene vs. scale');
    figure('Name', figureName);
    hold all;
    xlabel('$s$', 'Interpreter', 'LaTex'); 
    ylabel('$2 \Delta $', 'Interpreter', 'LaTex');

    plot(scale, slope, '-b');
    plot(scale, slope - confint99_lower, '-b');
    plot(scale, slope + confint99_upper, '-b');

    %errorbar( scale, slope, confint99_lower, confint99_upper, 'x', 'MarkerEdgeColor', 'k');

    %exactValue = -2;
    exactValue = -4;
    %exactValue = -0.1668575;        
    line([scale(1) ; scale(end)],[exactValue ; exactValue], 'Color','g', 'LineWidth', 1);

    hold off;

end



  