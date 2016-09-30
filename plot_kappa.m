cc = 0.9:0.0001:1.2;


kappa = zeros(1, numel(cc));


for kk=1:numel(cc)
   kappa(kk) = 6/(cc(kk)*(sqrt(12/cc(kk)) + 1)); 
    
end


figure; plot(kappa, cc);



S_test = zeros(1,numel(cc));

for kk=1:numel(cc)

    S_test(kk) = 1/(sqrt(12/cc(kk)) +1);
    
end

figure; plot(S_test, cc);