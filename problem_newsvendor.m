n = 1; % dimension of x 
Ain_x = [ 1; -1 ]; % feasible set of x 
bin_x = [ 100; 0 ]; % feasible set of x 

c = 1; % cost 
q = 5; % price 
n_2 = 2; % dimension of z
a = [ 1; -1]; 
B = [ 1 -1; 1 -1]; 
C = [ c-q; c ]; 
E = [ 0; -q ];
b = zeros(2,1);
r = size(B,1); 

% true distribution of xi 
m = 1; % dimension of xi
I_m = eye(m); 
xi_l = 0;
xi_u = 100;
Ain_xi = [ 1; -1 ]; % A matrix of support Xi
bin_xi = [ xi_u; -xi_l ]; % b matrix of support Xi 
dist_xi = makedist('Normal', 'mu', 20, 'sigma', 30); % mu = 20, 80
dist_xi = truncate(dist_xi, xi_l, xi_u); % true distribution

% PWA representation of f(x,xi) = max{f_1, .., f_K}
% f_k = xi'*A_k*x + b_k'*xi + c'_k*x + d_k
K = 2; 
A_k_matrix = zeros(K*m,n);
b_k_matrix = [ 0; -q ]; 
b_k_trans_matrix = [ 0; -q ]; 
b_k_cell = cell(2,1);
b_k_cell{1} = 0;
b_k_cell{2} = -q;
c_k_matrix = [ c - q; c ];
c_k_trans_matrix = [ c - q; c ]; 
c_k_cell = cell(2,1);
c_k_cell{1} = c - q;
c_k_cell{2} = c;
d_k_matrix = [0; 0];

% approx. the distribution of xi
N_true = 1e5; 
xi_vector_true = random(dist_xi, N_true, 1); 
xi_matrix_true = reshape(xi_vector_true, m, N_true)'; 

% true expected cost 
N_true_base = 250; 
N_true_set = N_true/N_true_base; 
Ain = [
    Ain_x zeros(size(Ain_x,1),N_true_base);
    kron(ones(N_true_base,1),c_k_trans_matrix) -kron(eye(N_true_base),ones(K,1));
    ];
bin = cell(N_true_set, 1);
for i_set = 1 : N_true_set 
    bin{i_set} = [
        bin_x;
        -kron(xi_matrix_true( N_true_base*(i_set-1) + 1 : N_true_base*i_set,: ), b_k_matrix) - kron(ones(N_true_base,1),d_k_matrix)
        ];
end
rho = 0.1; 
sol_new = [ kron(ones(1,N_true_set), zeros(n,1)) ];
y_new = zeros(n, N_true_set);
z_new = mean(sol_new + (1/rho) * y_new, 2);
fval_test = zeros(N_true_set,1);
err_r = 1;
err_s = 1;
H = rho * blkdiag( eye(n) , zeros(N_true_base) ); 
while (err_r > 1e-5 || err_s > 1e-5)
    for i = 1 : N_true_set
        sol_old = sol_new;
        f = [ 
            -rho * z_new + y_new(:,i); 
            1/N_true_base * ones(N_true_base,1) 
            ]; 
        sol = cplexqp(H,f,Ain,bin{i});
        sol_new(:,i) = sol(1:n);
        fval_test(i,1) = [ zeros(n,1); 1/N_true_base * ones(N_true_base,1) ]' * sol;
    end    
    z_old = z_new;
    z_new = mean(sol_new + (1/rho) * y_new, 2);   
    y_old = y_new;
    y_new = y_old + rho * (sol_new - z_new);
    r = sol_new - kron(ones(1,N_true_set), z_new);
    s = -rho * kron(kron(N_true_set,1),eye(n)) * ( z_new - z_old );
    err_r = norm(reshape(r,[],1))^2;
    err_s = norm(reshape(s,[],1))^2;
    
end
oos_true = mean(fval_test);

G = 1e3; % Big M
err_tol = 1e-5; % tolerance 
