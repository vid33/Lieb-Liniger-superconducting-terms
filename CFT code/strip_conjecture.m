clear;

uu = -5;
potential = 10;
interaction = 0;

kappa = 6/(sqrt(12)+1);

%uu=1;potential = -1; interaction=10;

fIn = sprintf('data_LLuu/lr_eigenvectors/uu=%dv=%dc=%deigenvalues_all.mat',uu, potential, interaction);

load(fIn, 'dim', 'eigvals');

entropy = zeros(1, numel(dim));
renyiEntropy = zeros(1, numel(dim)); alpha = 2;

%we go up to ""th smallest non-zero eigenvalue 
eigval_no = 10;

%where to start plots; lowest D = eigval_start + 3%
eigval_start =1;

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


figure;

plot(dim(eigval_start:end), -(1/6)*mu{3}.*(dim.^kappa) );
    