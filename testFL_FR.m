clear;

D=8;

uu = -5;
potential = 10;
interaction = 0;

%for calculating F and density matrices
bicg_maxiter=1000;
bicg_tolerance=1e-10;


    fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);
    load(fIn);

        FL_zero = zeros(D^2, 1);
        FR_zero = zeros(D^2, 1);
        
        [FL] = calculateF( R, Rkin, Q, matr, potential, mass, uu, interaction, FL_zero, D, bicg_maxiter, bicg_tolerance );
        
        
            [FR] = calculateF( R, Rkin, Q, matr, potential, mass, uu, interaction, FR_zero, D, bicg_maxiter, bicg_tolerance );

            
            matFL = transpose(reshape(FL, D, D));
            matFR = transpose(reshape(FR, D, D));
            
            HHtest = kron(transpose(matFL), inv(matr)*matFR) + kron(transpose(matFL), eye(D)) + kron(eye(D), inv(matr)*matFR);
            
            eigvals_test = sort(eig(HHtest));
            
            eigvals_test2 = eigvals_test - eigvals_test(1);
            
            TT = kron(R, conj(R)) + kron(Q, eye(D)) + kron(eye(D), conj(Q));
            
            eigvalsTT = -sort(eig(TT));
            
   [excitation_energies, H, Heigvectors] = calculateExcitations(0, R, Q, l, r, matr, 0, D, mass, interaction, potential, uu, 0);
         
            