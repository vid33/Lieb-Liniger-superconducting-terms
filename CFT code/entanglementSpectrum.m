clear;

D=64;

potential = 10;
interaction = 0;
uu = -5;

fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);

load(fIn);


lambda_sq = svd(matr);

ES = -log(sqrt(lambda_sq));


figure; plot((ES-ES(1))/(ES(2)-ES(1)), 'x');
