clear;

uu = -5;
potential = 10;
interaction = 0;

interpolationStart = 21; %position in array below where we interpolate to find central charge
eigvalPos = 2; %which eigval are we looking at; 2nd largest is meant to evaluate the central charge

dim = [ 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, ...
    31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, ...
    59, 60, 61, 62, 63, 64];

%dim = [ 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,31, 32];

%dim = [ 4, 8, 16, 20, 24, 32, 64];

entropy = zeros(1, numel(dim));
renyiEntropy = zeros(1, numel(dim)); alpha = 2;

%we go up to ""th smallest non-zero eigenvalue 
eigval_no = 30;

%where to start plots; lowest D = eigval_start + 4%
dim_start = 22;

mu_plot = cell(1, eigval_no);
dim_plot = cell(1, eigval_no);

for kk=1: eigval_no
    fIn = sprintf('data_LLuu/lr_eigenvectors/uu=%dv=%dc=%deigenvector=%d.mat',uu,  potential, interaction , kk+1);
    load(fIn, 'mu', 'dim');
    mu_plot{kk} = mu;
    dim_plot{kk} = dim;   
end

%logmu_sum= zeros(1, numel(dim));
%for kk=1: eigval_no
   
 %   logmu_sum = log(-1./mu_plot{4}) + log(-1./mu_plot{2});
   % logmu_sum = logmu_sum + log(-1./mu_plot{kk});
    
%end

figureName = sprintf('ln(corr) vs. ln(D)');
figure('Name', figureName);
        xlabel('$\log(D) $', 'Interpreter', 'LaTex'); 
        ylabel('$\log( \mu_I)$', 'Interpreter', 'LaTex');

eigval_end = eigval_no;  eigval_true = [0 1 0 1 0 1 0 1 0];
hold all;
for pp=2:eigval_end
    %if (eigval_true(pp) == 1)
    plot(log(dim_plot{pp}(dim_start:end)), log(-1./real(mu_plot{pp}(dim_start:end))), 'x');
    %end
end
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


