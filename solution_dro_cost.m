Ain_dro = [
    Ain_x zeros(size(Ain_x,1),1+N+size(Ain_xi,1)*N*K);
    blkdiag(dum_cell_x{:})*kron(ones(N,1),A_k_matrix)+kron(ones(N,1),c_k_trans_matrix) zeros(K*N,1) -kron(eye(N),ones(K,1)) blkdiag(dum_cell_gamma{:});
    kron(ones(N,1),A_k_matrix) -ones(m*N*K,1) zeros(m*N*K,N) -kron(eye(N*K),Ain_xi'); 
    -kron(ones(N,1),A_k_matrix) -ones(m*N*K,1) zeros(m*N*K,N) kron(eye(N*K),Ain_xi'); 
    ];
bin_dro = [
    bin_x
    -blkdiag(dum_cell_x{:})*kron(ones(N,1),b_k_matrix)-kron(ones(N,1),d_k_matrix);
    -kron(ones(N,1),b_k_matrix);
    kron(ones(N,1),b_k_matrix);
    ];
ub_dro = [
    inf * ones(n + 1 + N + size(Ain_xi,1)*N*K,1);
    ];
lb_dro = [ 
    -inf * ones(n,1); 
    0; 
    -inf * ones(N,1); 
    zeros(size(Ain_xi,1)*N*K,1);
    ];
f_dro = [ zeros(n,1); Was_dist; 1/N * ones(N,1); zeros(size(Ain_xi,1)*N*K,1) ];

sol = cplexlp(f_dro,Ain_dro,bin_dro,[],[],lb_dro,ub_dro);    

x = sol(1:n);
