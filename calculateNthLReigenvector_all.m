clear;

uu = -5;
potential = 10;
interaction = 0;

eigenvectorNo_max = 500; %store only the first ... eigenvectors

fOut = sprintf('data_LLuu/lr_eigenvectors/%dEIGVECSuu=%dv=%dc=%d.mat',eigenvectorNo_max, uu,  potential, interaction);


%dim = [ 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, ...
%    31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64];
%dim = 4:64;

dim = 2:64;

%dim = [ 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24];

%dim = [ 4, 8, 16, 20, 24, 32, 64];
mu = cell(eigenvectorNo_max, numel(dim));
%mu_tmp = zeros(1, numel(dim));

l2 = cell(eigenvectorNo_max, numel(dim));
r2 = cell(eigenvectorNo_max, numel(dim));

matl2 = cell(eigenvectorNo_max, numel(dim));
matr2 = cell(eigenvectorNo_max, numel(dim));


for k=1:numel(dim)
    
    fprintf('%d ', dim(k));
    if (rem(k,10)==0)
            fprintf('\n');
    end 
    
    fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', dim(k), dim(k), uu, potential, interaction);
 
    load(fIn, 'R', 'Q', 'D');
    
    eigenvectorNo_max_tmp = min(D^2, eigenvectorNo_max);
    
    TT = kron(Q, eye(D)) + kron(eye(D), conj(Q)) + kron(R, conj(R)); 
    
    [Vr, Dr] = eig(TT);
    [Vl, Dl] = eig(TT.'); Vl = conj(Vl);
    
    [Dr, Dr_index] = sort(diag(Dr));
    [Dl, Dl_index] = sort(diag(Dl));
    
    for zz=1:eigenvectorNo_max_tmp
     
        r2{zz, k} = Vr(:, Dr_index(zz));
        matr2{zz, k} = transpose(reshape(r2{zz, k}, D, D));
        
        l2{zz, k} = (Vl(:, Dl_index(zz)) )';
        matl2{zz, k} = transpose(reshape(l2{zz, k}, D, D));
        
        mu{zz, k} = Dl(zz);
        
    end  
    
 %   [ l2{k}, r2{k}, matl2{k}, matr2{k}, mu(k), r2_eigval_test] = lrEigenvector_subleading( R, Q, D, eigenvector);
    
    
end

clearvars TT Dr Dr_index Dl Dl_index Q R D;

fprintf('\n');





