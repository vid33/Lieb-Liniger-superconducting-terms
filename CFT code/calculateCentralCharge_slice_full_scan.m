%%PLOT log(\mu_I) vs. log(D) for any desired range of I Calculate entropies

clear;

potential = 10;
interaction = 0;
uu = -5;

symGauge = false;

scale = 0.1;

plot_delta_vs_scale99 = true;
plot_delta_vs_scale95 = true;

fIn = sprintf('data_LLuu/lr_eigenvectors/uu=%dv=%dc=%deigenvalues_all.mat', uu, potential, interaction);

load(fIn, 'dim', 'eigvals');

entropy = zeros(1, numel(dim));

centralCharge = 1; %expected central charge...

kappa = 6/( centralCharge*(sqrt(12/centralCharge)+1) );

%interpolate using nth eigalue 
interpol_eigval = 2;

%where to start plots; lowest D = eigval_start + 3%
eigval_start = 22;

%do interpolation starting for all ranges starting from eigval_start:dim(end) to eigval_start+scan_end:dim(end) 
scanDelta = 10;


mu = zeros(1, numel(dim));

%eigvals is sorted w.r.t. abs, uncomment this to sort using real part
for pp=1:numel(dim)
   eigvals{pp}   = -1*sort(-1*real(eigvals{pp}));
end

for pp=1:numel(dim)
    mu(pp) = eigvals{pp}(interpol_eigval);
end


for pp=1:numel(dim)
    
      fprintf('%d ', dim(pp));
        if (rem(dim(pp),10)==0)
                fprintf('\n');
        end 
    
    fIn =sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', dim(pp), dim(pp), uu, potential, interaction);
    load(fIn);
    if symGauge == true
        symmetricGauge;
    end
    
    D = dim(pp);
    
    TT = kron(R, conj(R)) + kron(Q, eye(D)) + kron(eye(D), conj(Q));
    
    expTTtmp = expm((-scale/mu(pp))*TT);
    
    clearvars TT;
    
    expTTtmp = reshape(transpose(expTTtmp), D, D, D, D);
    
    expTTtmp2 = expTTtmp;
    for kk=1:D
        for mm = 1:D
            expTTtmp2(:, kk, mm, :) = expTTtmp(:, mm,  kk, :);
        end
    end
    expTTtmp = expTTtmp2;
    
    clearvars expTTtmp2;
    
    expTT = reshape((expTTtmp), D^2, D^2);
    
    clearvars expTTtmp;
    
    %%%%%%%%
        D = dim(pp);
    
    fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);
    load(fIn); matl = eye(D);
   
    if symGauge == true
        symmetricGauge;
    end

    tmpmat = kron((transpose(matl^(0.5))), (transpose(matr))^(0.5))*(transpose(expTT))*kron((matl)^0.5, transpose(matr)^0.5);
    schmidt = svd(tmpmat);
    
    clearvars expTT tmpmat;

    entropy(pp) = -sum(schmidt.*log(schmidt));
    
end


fprintf('\n');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%ln(D)/ln(mu) vs. entropy

fitMeasure_entropyD = zeros(1, scanDelta);
slope_entropyD = zeros(1, scanDelta);
intercept_entropyD = zeros(1, scanDelta);
sse_entropyD = zeros(1, scanDelta);
rmse_entropyD = zeros(1, scanDelta);
rsquare_entropyD = zeros(1, scanDelta);
adjrsquare_entropyD = zeros(1, scanDelta);
confint95_entropyD = cell(1, scanDelta);
confint99_entropyD = cell(1, scanDelta);

confint95_entropyD_upper = zeros(1, scanDelta);
confint95_entropyD_lower = zeros(1, scanDelta);
confint99_entropyD_upper = zeros(1, scanDelta);
confint99_entropyD_lower = zeros(1, scanDelta);

fitMeasure_entropyMU = zeros(1, scanDelta);
slope_entropyMU = zeros(1, scanDelta);
intercept_entropyMU = zeros(1, scanDelta);
sse_entropyMU = zeros(1, scanDelta);
rmse_entropyMU = zeros(1, scanDelta);
rsquare_entropyMU = zeros(1, scanDelta);
adjrsquare_entropyMU = zeros(1, scanDelta);
confint95_entropyMU = cell(1, scanDelta);
confint99_entropyMU = cell(1, scanDelta);

confint95_entropyMU_upper = zeros(1, scanDelta);
confint95_entropyMU_lower = zeros(1, scanDelta);
confint99_entropyMU_upper = zeros(1, scanDelta);
confint99_entropyMU_lower = zeros(1, scanDelta);


for zz=1:scanDelta
    xxD = log(dim(eigval_start+zz-1:end));
    xxMU = real(log(-1./mu(eigval_start+zz-1:end)));
    yy = entropy(eigval_start+zz-1:end);
    %plot(xx, yy, '-');

    [fitobject, gof] = fit(transpose(xxD), transpose(yy), 'poly1');
    slope_entropyD(zz) = fitobject.p1;
    intercept_entropyD(zz) = fitobject.p2;
    sse_entropyD(zz) = gof.sse;
    rmse_entropyD(zz) = gof.rmse;
    rsquare_entropyD(zz) = gof.rsquare;
    adjrsquare_entropyD(zz) = gof.adjrsquare;
    
    confint95_entropyD{zz} = confint(fitobject, 0.95);
    confint99_entropyD{zz} = confint(fitobject, 0.9973);
    
    confint95_entropyD_lower(zz) = confint95_entropyD{zz}(1,1);
    confint95_entropyD_upper(zz) = confint95_entropyD{zz}(2,1);

    confint95_entropyD_lower(zz) = slope_entropyD(zz) - confint95_entropyD_lower(zz);
    confint95_entropyD_upper(zz) = confint95_entropyD_upper(zz) - slope_entropyD(zz);

    confint99_entropyD_lower(zz) = confint99_entropyD{zz}(1,1);
    confint99_entropyD_upper(zz) = confint99_entropyD{zz}(2,1);

    confint99_entropyD_lower(zz) = slope_entropyD(zz) - confint99_entropyD_lower(zz);
    confint99_entropyD_upper(zz) = confint99_entropyD_upper(zz) - slope_entropyD(zz);
    
    [fitobject, gof] = fit(transpose(xxMU), transpose(yy), 'poly1');
    slope_entropyMU(zz) = fitobject.p1;
    intercept_entropyMU(zz) = fitobject.p2;
    sse_entropyMU(zz) = gof.sse;
    rmse_entropyMU(zz) = gof.rmse;
    rsquare_entropyMU(zz) = gof.rsquare;
    adjrsquare_entropyMU(zz) = gof.adjrsquare;
    
    confint95_entropyMU{zz} = confint(fitobject, 0.95);
    confint99_entropyMU{zz} = confint(fitobject, 0.9973);
    
    confint95_entropyMU_lower(zz) = confint95_entropyMU{zz}(1,1);
    confint95_entropyMU_upper(zz) = confint95_entropyMU{zz}(2,1);

    confint95_entropyMU_lower(zz) = slope_entropyMU(zz) - confint95_entropyMU_lower(zz);
    confint95_entropyMU_upper(zz) = confint95_entropyMU_upper(zz) - slope_entropyMU(zz);

    confint99_entropyMU_lower(zz) = confint99_entropyMU{zz}(1,1);
    confint99_entropyMU_upper(zz) = confint99_entropyMU{zz}(2,1);

    confint99_entropyMU_lower(zz) = slope_entropyMU(zz) - confint99_entropyMU_lower(zz);
    confint99_entropyMU_upper(zz) = confint99_entropyMU_upper(zz) - slope_entropyMU(zz);

end


figureName = sprintf('ln(D) vs. entropy slopes');
figure('Name', figureName);
        xlabel('Initial D'); 
        ylabel('Slope and errors');
        
hold all;

if (plot_delta_vs_scale95 == true)
    plot( dim(eigval_start:eigval_start+scanDelta-1), slope_entropyD, '-b');
    plot( dim(eigval_start:eigval_start+scanDelta-1), slope_entropyD - confint95_entropyD_lower, '-b');
    plot( dim(eigval_start:eigval_start+scanDelta-1), slope_entropyD + confint95_entropyD_upper, '-b');
end

if (plot_delta_vs_scale99 == true)   
    plot( dim(eigval_start:eigval_start+scanDelta-1), slope_entropyD, '-b');
    plot( dim(eigval_start:eigval_start+scanDelta-1), slope_entropyD - confint99_entropyD_lower, '-r');
    plot( dim(eigval_start:eigval_start+scanDelta-1), slope_entropyD + confint99_entropyD_upper, '-r'); 
end

hold off;

exactValue = kappa/6;       
line([dim(eigval_start) ; dim(end) ],[exactValue ; exactValue], 'Color','g', 'LineWidth', 1);


figureName = sprintf('ln(mu) vs. entropy slopes');
figure('Name', figureName);
        xlabel('Initial D'); 
        ylabel('Slope and errors');
        
hold all;

if (plot_delta_vs_scale95 == true)
    plot( dim(eigval_start:eigval_start+scanDelta-1), slope_entropyMU, '-b');
    plot( dim(eigval_start:eigval_start+scanDelta-1), slope_entropyMU - confint95_entropyMU_lower, '-b');
    plot( dim(eigval_start:eigval_start+scanDelta-1), slope_entropyMU + confint95_entropyMU_upper, '-b');
end

if (plot_delta_vs_scale99 == true)   
    plot( dim(eigval_start:eigval_start+scanDelta-1), slope_entropyMU, '-b');
    plot( dim(eigval_start:eigval_start+scanDelta-1), slope_entropyMU - confint99_entropyMU_lower, '-r');
    plot( dim(eigval_start:eigval_start+scanDelta-1), slope_entropyMU + confint99_entropyMU_upper, '-r'); 
end

exactValue = centralCharge/3;       
line([dim(eigval_start) ; dim(end) ],[exactValue ; exactValue], 'Color','g', 'LineWidth', 1);

hold off;







%fitOut = polyfit(log(dim(interpolationStart:end)), log(1./mu(interpolationStart:end)), 1);

%fprintf('\nThe slope is %d, and the intercept %d\n', fitOut(1), fitOut(2));

%syms cCharge_symb;
%centralCharge = eval(solve( 6/(cCharge_symb*(sqrt(12/cCharge_symb)+1)) - abs(fitOut(1)), 'cCharge_symb' ) );

%fprintf('The central charge is %d\n', centralCharge);

%plotDelta = 0.05;
%xVal = (log(dim(1))-0.1):plotDelta:(log(dim(end))+0.1);
%yVal = fitOut(1).*xVal + fitOut(2);
%yVal = 1.331.*xVal - 0.58158;
%hold all;
%plot(xVal, yVal);
%hold off;

%figureName = sprintf('Entropy vs. ln(D)');
%figure('Name', figureName); 
%plot(log(dim), entropy, 'x');

%figureName = sprintf('Renyi-%d entropy vs. ln(D)', alpha);
%figure('Name', figureName); 
%plot(log(dim), renyiEntropy, 'x');


