%%PLOT log(\mu_I) vs. log(D) for any desired range of I. Calculate entropy,
%%plot S vs. log(mu_I), obtain central charge

clear;

uu = -5;
potential = 10;
interaction = 0;

%uu=1;potential = -1; interaction=10;

centralCharge = 1; %expected central charge...

kappa = 6/( centralCharge*(sqrt(12/centralCharge)+1) );

fIn = sprintf('data_LLuu/lr_eigenvectors/uu=%dv=%dc=%deigenvalues_all.mat',uu, potential, interaction);

load(fIn, 'dim', 'eigvals');

entropy = zeros(1, numel(dim));
renyiEntropy = zeros(1, numel(dim)); alpha = 2;

%we go up to ""th smallest non-zero eigenvalue 
eigval_no = 15;

%where to start plots; lowest D = eigval_start + 3%
eigval_start =3;

%end befor maximum available dimension, dim(end-offset)
offset = 0;

for k=1:numel(dim)  
    fprintf('%d ', dim(k));
    if (rem(k,10)==0)
            fprintf('\n');
    end 
    
    fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', dim(k), dim(k), uu, potential, interaction);
 
    load(fIn, 'matr');
    
    schmidt = svd(matr);
    entropy(k) = sum(-schmidt.*log(schmidt));
    
    renyiEntropy(k) = (1/(1-alpha))*log(sum(schmidt.^alpha));
     
    clearvars matr;    
end
fprintf('\n');


mu = cell(1, eigval_no);
for pp=1:eigval_no
    mu{pp} = zeros(1, numel(dim));
end

%eigvals is sorted w.r.t. abs, uncomment this to sort using real part
for pp=1:numel(dim)
   eigvals{pp}   = -1*sort(-1*real(eigvals{pp}));
end

for kk=1:eigval_no
    for pp=1:numel(dim)
        if(kk < dim(pp)^2)
            mu{kk}(pp) = eigvals{pp}(kk);
        end
    end
end

figureName = sprintf('ln(corr) vs. ln(D)');
figure('Name', figureName);
        xlabel('$\log(D) $', 'Interpreter', 'LaTex'); 
        ylabel('$\log( \mu_I)$', 'Interpreter', 'LaTex');
hold all;

fitMeasure = zeros(1, numel(eigval_no)-1);
slope = zeros(1, numel(eigval_no)-1);
intercept = zeros(1, numel(eigval_no)-1);
sse = zeros(1, numel(eigval_no)-1);
rmse = zeros(1, numel(eigval_no)-1);
rsquare = zeros(1, numel(eigval_no)-1);
adjrsquare = zeros(1, numel(eigval_no)-1 );
confint95 = zeros(1, numel(eigval_no)-1);
confint99 = zeros(1, numel(eigval_no)-1 );
    
for pp=2:(eigval_no)
    %if (pp==2 || pp==5 || pp==11 ) 
   % if (pp==5) 
    
    xx = real(log(dim(eigval_start:end-offset)));
    yy = real(log(-1./mu{pp}(eigval_start:end-offset)));
    plot(xx, yy, '-x');
    
    %  end
    [fitobject, gof] = fit(transpose(xx), transpose(yy), 'poly1');

    slope(pp-1) = fitobject.p1;
    intercept(pp-1) = fitobject.p2;
    sse(pp-1) = gof.sse;
    rmse(pp-1) = gof.rmse;
    rsquare(pp-1) = gof.rsquare;
    adjrsquare(pp-1) = gof.adjrsquare;
    
    confint95{pp-1} = confint(fitobject, 0.95);
    confint99{pp-1} = confint(fitobject, 0.99);
  
 end
%end
hold off;


[bestSSE, bestSSEIndex] = min(sse);
[bestRsquare, bestRsquareIndex] = max(rsquare);

%bestCorrLength = corrLength(bestSSEIndex);

bestSlope = slope(bestSSEIndex);
bestIntercept = intercept(bestSSEIndex);
bestConfInt95 = confint95{bestSSEIndex};
bestConfInt99 = confint99{bestSSEIndex};


%%%plot entropy
figureName = sprintf('entropy vs. ln(corr)');
figure('Name', figureName);
        xlabel('$\log( \mu_I) $', 'Interpreter', 'LaTex'); 
        ylabel('Entropy');
        
hold all;

fitMeasure_entropy = zeros(1, numel(eigval_no)-1);
slope_entropy = zeros(1, numel(eigval_no)-1 );
intercept_entropy = zeros(1, numel(eigval_no)-1 );
sse_entropy = zeros(1, numel(eigval_no)-1 );
rmse_entropy = zeros(1, numel(eigval_no)-1 );
rsquare_entropy = zeros(1, numel(eigval_no)-1 );
adjrsquare_entropy = zeros(1, numel(eigval_no)-1 );
confint95_entropy = zeros(1, numel(eigval_no)-1 );
confint99_entropy = zeros(1, numel(eigval_no)-1 );

for pp=2:(eigval_no )

    xx = real(log(-1./mu{pp}(eigval_start:end-offset) ));
    %xx = log(dim(eigval_start:end-offset));
    yy = entropy(eigval_start:end - offset);
    plot(xx, yy, 'x');

    [fitobject, gof] = fit(transpose(xx), transpose(yy), 'poly1');
    slope_entropy(pp-1) = fitobject.p1;
    intercept_entropy(pp-1) = fitobject.p2;
    sse_entropy(pp-1) = gof.sse;
    rmse_entropy(pp-1) = gof.rmse;
    rsquare_entropy(pp-1) = gof.rsquare;
    adjrsquare_entropy(pp-1) = gof.adjrsquare;
    
    confint95_entropy{pp-1} = confint(fitobject, 0.95);
    confint99_entropy{pp-1} = confint(fitobject, 0.99);

end

hold off;

%%%%%%%%%%%%%%plot entropy vs. ln(D)

figureName = sprintf('ln(D) vs. entropy');
figure('Name', figureName);
        xlabel('$\log( D) $', 'Interpreter', 'LaTex'); 
        ylabel('Entropy');
xx = log(dim(eigval_start:end));
yy = entropy(eigval_start:end);
plot(xx, yy, 'x');

[bestSSE_entropy, bestSSEIndex_entropy] = min(sse_entropy);
[bestRsquare_entropy, bestRsquareIndex_entropy] = max(rsquare_entropy);

bestSlope_entropy = slope_entropy(bestSSEIndex_entropy);
bestIntercept_entropy = intercept_entropy(bestSSEIndex_entropy);
bestConfInt95_entropy = confint95_entropy{bestSSEIndex_entropy};
bestConfInt99_entropy = confint99_entropy{bestSSEIndex_entropy};

%%% plot slope (\kappa) for all eigenvectors
confint99_upper = zeros(1, numel(eigval_no)-1);
confint99_lower = zeros(1, numel(eigval_no)-1);

for pp=2:eigval_no
    
    confint99_lower(pp-1) = confint99{pp-1}(1,1);
    confint99_upper(pp-1) = confint99{pp-1}(2,1);

end

confint99_lower = slope - confint99_lower;
confint99_upper = confint99_upper - slope;

eigenvalue_vec = 2:eigval_no;

%figure; plot(eigenvalue_vec, confint99_lower, 'x');

%hold all;



figureName = sprintf('kappa vs. eigenvalue');
figure('Name', figureName);
        xlabel('Eigenvalue');
        %ylabel('$\kappa $ (slope of $\log(\mu_I)$ vs. $\log(D)$)', 'Interpreter', 'LaTex'); 
        ylabel('$\kappa$ (slope of log(\mu_I) vs. log(D))'); 
%plot(eigenvalue_vec, confint99_upper, 'x');
hold('on');

plot(eigenvalue_vec, slope, 'k-');

errorbar( eigenvalue_vec, slope, confint99_lower, confint99_upper, 'b', 'Marker', 'none', 'LineStyle', 'none' );

exactValue = kappa;       
line([eigenvalue_vec(1) ; eigenvalue_vec(end) ],[exactValue ; exactValue], 'Color','g', 'LineWidth', 2);

hold off;

%see if we can say something about ratios

tmp_up = mu{5}(eigval_start:end-offset);
tmp_down = mu{2}(eigval_start:end-offset);

figure; plot(dim(eigval_start:end-offset), tmp_up./tmp_down);


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


figure; plot(dim(eigval_start:end-offset), mu{4}(eigval_start:end-offset)./mu{3}(eigval_start:end-offset), 'x-');


