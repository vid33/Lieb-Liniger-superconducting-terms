schmidt = zeros(D, N);

energyUpToX_FL = zeros(1, N);
energyUpToX_FR = zeros(1, N);
sum_FL_FR = zeros(1, N);
entropy = zeros(1, N);

%Rexpectation = zeros(1,N);


for kk=1:N
   
    XY = transpose(matl_cell{kk})^0.5*matr_cell{kk}*transpose(matl_cell{kk})^0.5;
    
    schmidt(:,kk) = svd(XY);
     energyUpToX_FL(kk) = trace(transpose(matFL{kk})*matr_cell{kk}); 
     energyUpToX_FR(kk) = trace(transpose(matl_cell{kk})*matFR{kk}); 
     sum_FL_FR(kk) = energyUpToX_FL(kk) + energyUpToX_FR(kk);
     %Rexpectation(kk) = trace(transpose(matl{kk})*R{kk}*matr{kk});
     entropy(kk) = sum(-schmidt(:, kk).*log(schmidt(:, kk)));
     
     fprintf('Schmidt sum is %d\n', sum(schmidt(:,kk)));
     
end

%Rexpectation = angle(Rexpectation);

figure; plot(position, transpose(schmidt));

figure; 
hold all;
plot(position, energyUpToX_FL);
plot(position, energyUpToX_FR);
plot(position, sum_FL_FR);
hold off;

%figure; plot(position(2:N-1), Rexpectation(2:N-1));
figure; plot(position, entropy);