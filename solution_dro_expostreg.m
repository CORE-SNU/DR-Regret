L = 0; 
UB = inf;
LB = -inf;

while UB - LB > err_tol 

    L = L + 1; 

    if L == 1 

        % initialization of w 
        z = zeros(N*n_2, 1); 
        q_plus = zeros(N*m, 1); 
        q_minus = zeros(N*m, 1); 

        % initialization of nu 
        nu = zeros(r*N, 1); 
    
        % initialization of constraints in MP
        Ain_master = [ 
            Ain_x zeros(size(Ain_x, 1), 2); 
            (1/N)*(kron(ones(1,N), eye(r))*nu)'*C -(1/N)*sum(q_plus + q_minus) -1; 
            ];
        bin_master = [
            bin_x;
            (1/N)*a'*kron(ones(1,N), eye(n_2))*z - (1/N)*(xi_hat_vec+q_plus-q_minus)'*kron(I_N, E')*nu - (1/N)*b'*kron(ones(1,N), eye(r))*nu
            ];
        lb_master = [
            zeros(n,1);
            0;
            -inf;
            ];

        % objective function of MP
        f_master = [ zeros(n,1); Was_dist; 1 ];    

    else

        % cut
        Ain_master = [
            Ain_master;
            (1/N)*(kron(ones(1,N), eye(r))*nu)'*C -(1/N)*sum(q_plus + q_minus) -1;             
            ];
        bin_master = [
            bin_master;
            (1/N)*a'*kron(ones(1,N), eye(n_2))*z-(1/N)*(xi_hat_vec+q_plus-q_minus)'*kron(I_N, E')*nu-(1/N)*b'*kron(ones(1,N), eye(r))*nu            
            ];

    end
    
    %% master problem     
    [ sol, LB ] = cplexlp(f_master,Ain_master,bin_master,[],[],lb_master);    
    x = sol(1 : n);
    lambda = sol(n + 1);
    eta = sol(end);

    %% subproblem
    Ain_sub = [
        zeros(N*size(Ain_x,1), r*N + n_2*N) kron(I_N,Ain_x) zeros(N*size(Ain_x,1), N*n_2 + 2*N*m + N*n_2 + N*r);
        zeros(N*size(Ain_xi,1), r*N + n_2*N + N*n + N*n_2) kron(I_N,Ain_xi) -kron(I_N,Ain_xi) zeros(N*size(Ain_xi,1), N*n_2 + N*r);
        zeros(N*size(C,1), r*N + n_2*N) kron(I_N,C) -kron(I_N,B) kron(I_N,E) -kron(I_N,E) zeros(N*size(C,1), N*n_2 + N*r);
        zeros(N*size(C,1), r*N + n_2*N + N*n + N*n_2) kron(I_N,E) -kron(I_N,E) -kron(I_N,B) zeros(N*size(C,1), N*r);
        zeros(N*size(B',1), r*N + n_2*N + N*n + N*n_2 + 2*N*m + N*n_2) kron(I_N,B');
        -G*eye(N*r) zeros(N*r,N*n_2 + N*n + N*n_2) -kron(I_N,E) kron(I_N,E) kron(I_N,B) zeros(N*r);
        G*eye(N*r) zeros(N*r,N*n_2 + N*n + N*n_2 + 2*N*m + N*n_2) eye(N*r);
        zeros(N*n_2) -G*eye(N*n_2) zeros(N*n_2, N*n + N*n_2 + 2*N*m + N*n_2) -kron(I_N,B');
        zeros(N*n_2) G*eye(N*n_2) zeros(N*n_2, N*n + N*n_2 + 2*N*m) eye(N*n_2) eye(N*n_2,N*r);
        ];
    bin_sub = [
        kron(ones(N,1), bin_x);
        kron(ones(N,1), bin_xi)-kron(I_N, Ain_xi)*xi_hat_vec;
        -kron(I_N,E)*xi_hat_vec-kron(ones(N,1), b);
        -kron(ones(N,1), C*x)-kron(I_N,E)*xi_hat_vec-kron(ones(N,1),b); 
        kron(ones(N,1),a);
        kron(ones(N,1), C*x)+kron(I_N,E)*xi_hat_vec+kron(ones(N,1),b); 
        G*ones(N*r,1); 
        -kron(ones(N,1),a);
        G*ones(N*n_2,1);         
        ];
    f_sub = [
        zeros(r*N + N*n_2 + N*n, 1);
        (1/N) * kron(ones(N,1),a);
        (1/N) * lambda * ones(N*m,1);
        (1/N) * lambda * ones(N*m,1);
        -(1/N) * kron(ones(N,1),a);
        zeros(r*N,1)
        ]; 
    lb_sub = [ 
        zeros(r*N + N*n_2 + N*n + N*n_2 + 2*N*m + N*n_2 + r*N, 1);
        ];
    ub_sub = [
        ones(r*N + N*n_2, 1);
        inf*ones(N*n + N*n_2 + 2*N*m + N*n_2 + r*N, 1);
        ];
    ctype_sub = [ repmat('B',1,r*N + N*n_2) repmat('C',1,N*n + N*n_2 + 2*N*m + N*n_2 + r*N) ];
    
    [ sol_sub, fval_sub ] = cplexmilp(f_sub,Ain_sub,bin_sub,[],[],[],[],[],lb_sub,ub_sub,ctype_sub);

    z = sol_sub(r*N + N*n_2 + N*n + 1 : r*N + N*n_2 + N*n + N*n_2);
    q_plus = sol_sub(r*N + N*n_2 + N*n + N*n_2 + 1 : r*N + N*n_2 + N*n + N*n_2 + N*m);
    q_minus = sol_sub(r*N + N*n_2 + N*n + N*n_2 + N*m + 1 : r*N + N*n_2 + N*n + N*n_2 + N*m + N*m);
    nu = sol_sub( r*N + N*n_2 + N*n + N*n_2 + 2*N*m + N*n_2 + 1 : r*N + N*n_2 + N*n + N*n_2 + 2*N*m + N*n_2 + r*N);

    UB = Was_dist * lambda - fval_sub;
       
end
