clear;

uu = -5;
potential = 10;
interaction = 0;

eigenvector = 9; %nth eigenvector

fOut = sprintf('data_LLuu/lr_eigenvectors/uu=%dv=%dc=%deigenvector=%d.mat',uu,  potential, interaction , eigenvector);


dim = [ 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, ...
    31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64];


%dim = [ 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24];

%dim = [ 4, 8, 16, 20, 24, 32, 64];
mu = zeros(1, numel(dim));
mu_tmp = zeros(1, numel(dim));

l2 = cell(1, numel(dim));
r2 = cell(1, numel(dim));

matl2 = cell(1, numel(dim));
matr2 = cell(1, numel(dim));


for k=1:numel(dim)
    
    fprintf('%d ', dim(k));
    if (rem(k,10)==0)
            fprintf('\n');
    end 
    
    fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', dim(k), dim(k), uu, potential, interaction);
 
    load(fIn, 'R', 'Q', 'D', 'matr', 'l', 'r', 'occupationDensity');
    
    [ l2{k}, r2{k}, matl2{k}, matr2{k}, mu(k), r2_eigval_test] = lrEigenvector_subleading( R, Q, D, eigenvector);
    
    
end

fprintf('\n');





