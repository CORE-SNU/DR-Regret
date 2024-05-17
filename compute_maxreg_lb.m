Ain_sub = [
        zeros(size(Ain_x,1), r*N + n_2*N) Ain_x zeros(size(Ain_x,1), N*n_2 + 2*N*m + N*n_2 + N*r);
        zeros(N*size(Ain_xi,1), r*N + n_2*N + n + N*n_2) kron(I_N,Ain_xi) -kron(I_N,Ain_xi) zeros(N*size(Ain_xi,1), N*n_2 + N*r);
        zeros(N*size(C,1), r*N + n_2*N) kron(ones(N,1),C) -kron(I_N,B) kron(I_N,E) -kron(I_N,E) zeros(N*size(C,1), N*n_2 + N*r);
        zeros(N*size(C,1), r*N + n_2*N + n + N*n_2) kron(I_N,E) -kron(I_N,E) -kron(I_N,B) zeros(N*size(C,1), N*r);
        zeros(N*size(B',1), r*N + n_2*N + n + N*n_2 + 2*N*m + N*n_2) kron(I_N,B');
        -G*eye(N*r) zeros(N*r,N*n_2 + n + N*n_2) -kron(I_N,E) kron(I_N,E) kron(I_N,B) zeros(N*r);
        G*eye(N*r) zeros(N*r,N*n_2 + n + N*n_2 + 2*N*m + N*n_2) eye(N*r);
        zeros(N*n_2) -G*eye(N*n_2) zeros(N*n_2, n + N*n_2 + 2*N*m + N*n_2) -kron(I_N,B');
        zeros(N*n_2) G*eye(N*n_2) zeros(N*n_2, n + N*n_2 + 2*N*m) eye(N*n_2) eye(N*n_2,N*r);
        zeros(1,r*N + n_2*N + n + N*n_2) (1/N)*ones(1,N*m) (1/N)*ones(1,N*m) zeros(1,N*n_2 + N*r);    
        ];
bin_sub = [
    bin_x;
    kron(ones(N,1), bin_xi)-kron(I_N, Ain_xi)*xi_hat_vec;
    -kron(I_N,E)*xi_hat_vec-kron(ones(N,1), b);
    -kron(ones(N,1), C*x)-kron(I_N,E)*xi_hat_vec-kron(ones(N,1),b); 
    kron(ones(N,1),a);
    kron(ones(N,1), C*x)+kron(I_N,E)*xi_hat_vec+kron(ones(N,1),b); 
    G*ones(N*r,1); 
    -kron(ones(N,1),a);
    G*ones(N*n_2,1);
    Was_dist;   
    ];
f_sub = [
    zeros(r*N + N*n_2 + n, 1);
    (1/N) * kron(ones(N,1),a);
    zeros(N*m,1);
    zeros(N*m,1);
    -(1/N) * kron(ones(N,1),a);
    zeros(r*N,1) 
    ]; 
lb_sub = [ 
    zeros(r*N + N*n_2 + n + N*n_2 + 2*N*m + N*n_2 + r*N, 1);
    ];
ub_sub = [
    ones(r*N + N*n_2, 1);
    inf*ones(n + N*n_2 + 2*N*m + N*n_2 + r*N, 1);
    ];
ctype_sub = [ repmat('B',1,r*N + N*n_2) repmat('C',1,n + N*n_2 + 2*N*m + N*n_2 + r*N) ];

[ ~, fval_sub ] = cplexmilp(f_sub,Ain_sub,bin_sub,[],[],[],[],[],lb_sub,ub_sub,ctype_sub);
lb_maxreg = - fval_sub;
       