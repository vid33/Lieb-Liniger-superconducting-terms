clear;

potential = 10;
interaction = 0;
uu = -5;
mass = 0.5;

dimStart = 4;  dimEnd = 20;

dim = dimStart:1:dimEnd;

centralCharge = 1; %expected central charge...

kappa = 6/( centralCharge*(sqrt(12/centralCharge)+1) );

scale  =1;

%l2 = cell(eigenvectorMax, numel(dim));
%r2 = cell(eigenvectorMax, numel(dim));

%matl2 = cell(eigenvectorMax, numel(dim));
%matr2 = cell(eigenvectorMax, numel(dim));
%mu = cell(eigenvectorMax, numel(dim));

expTT = cell(1, numel(dim));

entropyL = zeros(1,  numel(dim));
entropyR = zeros(1,  numel(dim));
entropyLR = zeros(1,  numel(dim));
entropy0 = zeros(1,  numel(dim));
mu2 = zeros(1, numel(dim));

for pp=1:numel(dim)
    
      fprintf('%d ', dim(pp));
        if (rem(dim(pp),10)==0)
                fprintf('\n');
        end 
    
        fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', dim(pp), dim(pp), uu, potential, interaction);
    %fIn = sprintf('data/D=%d/cpxD=%dv=%dc=%d.mat',dim(pp), dim(pp), -potential, interaction);
    load(fIn); matl = eye(D);
    
    D = dim(pp);
    
    fIn = sprintf('data_LLuu/D=%d/boundary/boundaryL_D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);
    load(fIn, 'VL', 'matVL', 'VLsqrt', 'norm');
    
    fIn = sprintf('data_LLuu/D=%d/boundary/boundaryR_D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);
    load(fIn, 'VR', 'matVR', 'VRsqrt');
    

    TT = kron(R, conj(R)) + kron(Q, eye(D)) + kron(eye(D), conj(Q));
    eigTT = sort(eig(TT));
    mu2(pp) = eigTT(2);
    
    vl = VL*expm((-scale/mu2(pp))*TT);
    vr = expm((-scale/mu2(pp))*TT)*VR;
    matvl = transpose(reshape(vl, D, D));
    matvr = transpose(reshape(vr, D, D));
   
    schmidt0 = svd(matr);
    entropy0(pp) = sum(-schmidt0.*log(schmidt0));
    
    %schmidtL = svd(transpose(matvl)*matr);
     schmidtL = svd(transpose(matvl)^0.5*matr*transpose(matvl)^0.5);
    entropyL(pp) = -sum(schmidtL.*log(schmidtL));
    schmidtR = svd(matvr);
    entropyR(pp) = -sum(schmidtR.*log(schmidtR));
    
    %norm_tmp = vl*vr;
    %vl = vl/norm_tmp;
    %matvl =matvl/norm_tmp;
    
    %schmidtLR = svd(transpose(matvl)^0.5*matvr*transpose(matvl)^0.5);
    %entropyLR(pp) = -sum(schmidtLR.*log(schmidtLR));

end

fprintf('\n');

figureName0 = sprintf('Entropy at infty');
figureNameL = sprintf('Entropy L');
figureNameR = sprintf('Entropy R');

figure('Name', figureName0);
plot(log(1./mu2), entropy0, 'x');

figure('Name', figureNameL);
plot(log(-1./mu2), entropyL, 'x');

figure('Name', figureNameR);
plot(log(-1./mu2), entropyR, 'x');
%figure; plot(log(-1./mu2), entropyLR, 'x');

