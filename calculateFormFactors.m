clear; ff_tic = tic;

D=16;
%potential =  -1;
%interaction = 10;
%uu = 1;

uu = -5;
potential = 10;
interaction = 0;


DSF = false;
%in case we want to store all the eigenvectors. Takes up a lot of memory for high D
store = false;

fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);
load(fIn);
occupationDensity = real(occupationDensity);

p =0*occupationDensity;

if (DSF == true)
    fOut = sprintf('data_LLuu/D=%d/DSF_form_factors/DSF_form_factorsD=%duu=%dv=%dc=%dp=%d.mat', ...
    D , D, uu, potential, interaction , p);
elseif (DSF == false)
    fOut = sprintf('data_LLuu/D=%d/form_factors/form_factorsD=%duu=%dv=%dc=%dp=%d.mat', ...
        D , D, uu, potential, interaction , p);
end

fprintf('Calculating Tminus...\n');
Tminus_tic = tic;
Tminus = -kron(Q, eye(D))- kron( eye(D), conj(Q)) - kron(R, conj(R)) - kron(r,l) - 1i*p*eye(D^2);
%Tzero = -T - kron(r,l);
%Tminus =  Tzero + 1i*p*eye(D^2, D^2);
Tminus_duration = toc(Tminus_tic);
fprintf('...took %d sec\n', Tminus_duration);

%[excitation_energies, H, Heigvectors] = calculateExcitations(p, R, Q, l, r, matr, Tminus_in, D, mass, interaction, potential, uu, FOR_PLOT)
[excitation_energies, H, Heigvectors] = calculateExcitations(p, R, Q, l, r, matr, Tminus, D, mass, interaction, potential, uu, 0);

C = cell(1,D^2);
if (DSF == false)
    C1P = cell(1, D^2); % C 1-particle, hole contribution to green, C itself is the particle contribution
end

Tminus = Tminus +2i*p*eye(D^2, D^2); %actually using Tplus here, it's huge for large D so don't wanna copy it

if (DSF == true)
    tmp_left = reshape(R'*R, 1, D^2)/Tminus;
elseif (DSF == false)
    tmp_left=reshape(R', 1, D^2)/Tminus;
    tmp_left_1P = reshape(R, 1, D^2)/Tminus; %hole contribution
    tmp_left_1P = reshape(tmp_left_1P, D, D);
end
tmp_left = reshape(tmp_left, D, D); %'missing' transpose not a mistake! 

clearvars Tminus;

matr_sqrt = sqrt(matr);
matrinv_sqrt = matr^(-1/2);


if (store == true)
    W = cell(1, D^2);
    V = cell(1, D^2);
end

 for k = 1: D^2
        fprintf('%d ', k);
        if (rem(k,10)==0)
            fprintf('\n');
        end  
        if (store == true)
            W{k} = transpose(reshape(Heigvectors(:, k), D, D))*matrinv_sqrt;
            V{k} = -(R')*W{k};
        elseif (store == false)
            W = transpose(reshape(Heigvectors(:, k), D, D))*matrinv_sqrt;
            V = -(R')*W;
        end
        
        if (store == true)
            if (DSF == true)
                C{k} = trace(R*matr*W{k}') + trace(tmp_left*matr*V{k}') + trace(tmp_left*R*matr*W{k}');
            elseif (DSF == false)
                C{k} = trace(matr*W{k}') + trace(tmp_left*matr*V{k}') + trace(tmp_left*R*matr*W{k}');
                C1P{k} = trace(tmp_left_1P*matr*V{k}') + trace(tmp_left_1P*R*matr*W{k}'); %hole
            end
        end
        
        if (store == false)
            if (DSF == true)
                C{k} = trace(R*matr*W') + trace(tmp_left*matr*V') + trace(tmp_left*R*matr*W');
            elseif (DSF == false)
                C{k} = trace(matr*W') + trace(tmp_left*matr*V') + trace(tmp_left*R*matr*W');
                C1P{k} = trace(tmp_left_1P*matr*V') + trace(tmp_left_1P*R*matr*W'); %hole
            end
        end
        
 end
 fprintf('\n');

 %save disk space for large D
 if (store == true)
    clearvars H Heigvectors tmp_left tmp_left_1P matrinv_sqrt; 
 elseif (store == false)  
    clearvars H Heigvectors tmp_left tmp_left_1P W V matrinv_sqrt;
 end
 
 ff_duration = toc(ff_tic);
 fprintf('calculateFormFactors took %f sec\n', ff_duration);
 
 %save(fOut);

