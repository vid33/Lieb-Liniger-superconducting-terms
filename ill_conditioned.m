illA=randmat(4, 1e60);

illA_I = illA/illA;


illA_inv = invillco(illA);

illA_I2 = prodKL(illA_inv, illA, 10, 1);

A = rand(D);

A_I = A/A;

A_inv = invillco(A);

A_I2 = prodKL(A, A_inv, 3, 1);