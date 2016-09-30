function M = mass_matrix_ode( t, matr, D)

    M = sparse(transpose(kron(eye(D), matr)));


end

