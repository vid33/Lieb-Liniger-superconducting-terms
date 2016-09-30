clear;

uu = -5;
potential = 10;
interaction = 0;

eigenvectorNo_max = 15;  
fIn = sprintf('data_LLuu/lr_eigenvectors/%dEIGVECSuu=%dv=%dc=%d.mat',eigenvectorNo_max, uu,  potential, interaction);
load(fIn, 'dim', 'mu', 'matl2', 'matr2');



currentDimension = 4;


entropy = zeros(1, eigenvectorNo_max);
XY  = cell(1, eigenvectorNo_max);
schmidt = cell(1, eigenvectorNo_max);
schmidt_sum = zeros(1, eigenvectorNo_max);
schmidt_min = zeros(1, eigenvectorNo_max);

for kk=1:eigenvectorNo_max
   
    %XY = transpose(matl2{kk, currentDimension-1})^0.5*matr2{kk, currentDimension - 1}^0.5;
    matl2_t = transpose(matl2{kk, currentDimension-1});
    XY = matl2_t^0.5*matr2{kk, currentDimension - 1}*matl2_t^0.5;
    
    schmidt{kk} = svd(XY);
    schmidt_min(kk) = min(schmidt{kk});
    schmidt_sum(kk) = sum(schmidt{kk});
    schmidt{kk} = schmidt{kk}/sum(schmidt{kk});
    entropy(kk) = sum(-schmidt{kk}.*log(schmidt{kk}));
    
    
end

figure; plot(1:eigenvectorNo_max, schmidt_min,'x-');
%figure; plot(1:eigenvectorNo_max, entropy);