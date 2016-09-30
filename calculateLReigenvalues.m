clear;

uu = 1;
potential = -1;
interaction = 10;

eigenvector = 1; %nth eigenvector

fOut = sprintf('data_LLuu/lr_eigenvectors/uu=%dv=%dc=%deigenvalues_all.mat', uu, potential, interaction);


%dim = [ 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, ...
%    31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64];

dim = [2, 3, 4, 5, 6, 7, 8, 9, 10 , 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 22];

%dim = [ 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24];

%dim = [ 4, 8, 16, 20, 24, 32, 64];

eigvals = cell(1, numel(dim));

for k=1:numel(dim)
    
    fprintf('%d ', dim(k));
    if (rem(k,10)==0)
            fprintf('\n');
    end 
    
    fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', dim(k), dim(k), uu, potential, interaction); 
 
    load(fIn, 'R', 'Q', 'D', 'matr', 'l', 'r', 'occupationDensity');
    
    TT = kron(R, conj(R)) + kron(Q, eye(D)) + kron(eye(D), conj(Q));
    
    eigvals{k}= sort(eig(TT));
    
    
end

clearvars('TT');

fprintf('\n');





