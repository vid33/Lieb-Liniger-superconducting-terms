function [ out] = cpxrand(D1, D2)



tmp_real = (rand(D1,D2)-0.5);
tmp_cpx = 1i*(rand(D1, D2) - 0.5);
%tmp_real = (rand(D1,D2));
%tmp_cpx = i*(rand(D1, D2));

out=tmp_real+tmp_cpx;



end

