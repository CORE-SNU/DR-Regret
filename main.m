clear 

%% Parameters 
problem_newsvendor % problem, distribution 
n_simu = 100; % simulation runs 
N = 10; % sample size 
I_N = eye(N); 
eps_set = norm(xi_u - xi_l, 1) * [ 0 1e-4 1e-3 1e-2 1e-1 1 ]; % Wasserstein ball's radii
n_eps = size(eps_set,2); 

%% Scen. gen.
sample_index = cell(n_simu,1);
for i_simu = 1 : n_simu
    sample_index{i_simu} = randi([1 N_true], N, 1); 
end

%% Sol. 
% dro_exantereg 
sol_dro_exantereg = cell(n_simu,n_eps); 
maxreg_ub_dro_exantereg = zeros(n_simu,n_eps); 
for i_simu = 1 : n_simu 
    xi_hat_matrix = xi_matrix_true(sample_index{i_simu},:); 
    xi_hat_vec = reshape(xi_hat_matrix',[],1);    
    for i_eps = 1 : n_eps
        Was_dist = eps_set(i_eps);
        solution_dro_exantereg
        sol_dro_exantereg{i_simu,i_eps} = x;
        maxreg_ub_dro_exantereg(i_simu,i_eps) = UB;
    end
end
maxreg_ub_dro_exantereg = mean(maxreg_ub_dro_exantereg);

% dro_expostreg
sol_dro_expostreg = cell(n_simu,n_eps);
for i_simu = 1 : n_simu 
    xi_hat_matrix = xi_matrix_true(sample_index{i_simu},:); 
    xi_hat_vec = reshape(xi_hat_matrix',[],1);    
    for i_eps = 1 : n_eps
        Was_dist = eps_set(i_eps);
        solution_dro_expostreg
        sol_dro_expostreg{i_simu,i_eps} = x;
    end
end

% dro_cost
sol_dro_cost = cell(n_simu,n_eps);
for i_simu = 1 : n_simu 
    xi_hat_matrix = xi_matrix_true(sample_index{i_simu},:); 
    xi_hat_vec = reshape(xi_hat_matrix',[],1); 
    dum_cell_x = cell(N,1);
    dum_cell_gamma = cell(N,1);
    for i = 1 : N 
        dum_cell_x{i} = kron(eye(K),xi_hat_matrix(i,:));
        dum_cell_gamma{i} = kron(eye(K),(bin_xi - Ain_xi * xi_hat_matrix(i,:)')');
    end
    for i_eps = 1 : n_eps
        Was_dist = eps_set(i_eps);
        solution_dro_cost
        sol_dro_cost{i_simu,i_eps} = x;
    end
end

%% UB on max. exantereg
% dro_expostreg
maxreg_ub_dro_expostreg = zeros(n_simu,n_eps);
for i_simu = 1 : n_simu
    xi_hat_matrix = xi_matrix_true(sample_index{i_simu},:); 
    xi_hat_vec = reshape(xi_hat_matrix',[],1);            
    for i_eps = 1 : n_eps
        x = sol_dro_expostreg{i_simu,i_eps}; 
        Was_dist = eps_set(i_eps);
        compute_maxreg_ub
        maxreg_ub_dro_expostreg(i_simu,i_eps) = ub_maxreg;       
    end
end
maxreg_ub_dro_expostreg = mean(maxreg_ub_dro_expostreg);        

% dro_cost
maxreg_ub_dro_cost = zeros(n_simu,n_eps);
for i_simu = 1 : n_simu
    xi_hat_matrix = xi_matrix_true(sample_index{i_simu},:); 
    xi_hat_vec = reshape(xi_hat_matrix',[],1);                
    for i_eps = 1 : n_eps
        x = sol_dro_cost{i_simu,i_eps};
        Was_dist = eps_set(i_eps); 
        compute_maxreg_ub
        maxreg_ub_dro_cost(i_simu,i_eps) = ub_maxreg;       
    end
end
maxreg_ub_dro_cost = mean(maxreg_ub_dro_cost); 

%% LB on max. exantereg
% dro_exantereg
maxreg_lb_dro_exantereg = zeros(n_simu,n_eps);        
for i_simu = 1 : n_simu    
    xi_hat_matrix = xi_matrix_true(sample_index{i_simu},:); 
    xi_hat_vec = reshape(xi_hat_matrix',[],1);    
    for i_eps = 1 : n_eps
        x = sol_dro_exantereg{i_simu,i_eps}; 
        Was_dist = eps_set(i_eps);
        compute_maxreg_lb
        maxreg_lb_dro_exantereg(i_simu,i_eps) = lb_maxreg;        
    end
end
maxreg_lb_dro_exantereg = mean(maxreg_lb_dro_exantereg);

% dro_expostreg 
maxreg_lb_dro_expostreg = zeros(n_simu,n_eps);
for i_simu = 1 : n_simu    
    xi_hat_matrix = xi_matrix_true(sample_index{i_simu},:); 
    xi_hat_vec = reshape(xi_hat_matrix',[],1);      
    for i_eps = 1 : n_eps
        x = sol_dro_expostreg{i_simu,i_eps}; 
        Was_dist = eps_set(i_eps);
        compute_maxreg_lb
        maxreg_lb_dro_expostreg(i_simu,i_eps) = lb_maxreg;
    end
end
maxreg_lb_dro_expostreg = mean(maxreg_lb_dro_expostreg);

% dro_cost 
maxreg_lb_dro_cost = zeros(n_simu,n_eps);
for i_simu = 1 : n_simu    
    xi_hat_matrix = xi_matrix_true(sample_index{i_simu},:); 
    xi_hat_vec = reshape(xi_hat_matrix',[],1);
    for i_eps = 1 : n_eps
        x = sol_dro_cost{i_simu,i_eps};
        Was_dist = eps_set(i_eps); 
        compute_maxreg_lb
        maxreg_lb_dro_cost(i_simu,i_eps) = lb_maxreg;
    end
end
maxreg_lb_dro_cost = mean(maxreg_lb_dro_cost);

%% OOS reg 
N_true_base = 250; 
N_true_set = N_true/N_true_base; 
f_oos = 1/N_true_base * ones(N_true_base,1); 
xi_vector_true_reshape = reshape(xi_vector_true, N_true_base*m, N_true_set);    
Ain_oos = -kron(eye(N_true_base),ones(K,1));

% dro_exantereg
fval_oos_dro_exantereg = zeros(N_true_set,n_eps);
oos_dro_exantereg = zeros(n_simu,n_eps);
for i_simu = 1 : n_simu         
    for i_set = 1 : N_true_set                 
        xi_hat_matrix_true = reshape( xi_vector_true_reshape(:,i_set), m, N_true_base )';
        dum_cell_x_true = cell(N_true_base, 1);
        for i_base = 1 : N_true_base
            dum_cell_x_true{i_base} = kron(eye(K), xi_hat_matrix_true(i_base,:)); 
        end    
        for i_eps = 1 : n_eps 
            x_oos = sol_dro_exantereg{i_simu,i_eps};        
            bin_oos = -blkdiag(dum_cell_x_true{:})*kron(ones(N_true_base,1),b_k_matrix) - kron(ones(N_true_base,1), c_k_trans_matrix*x_oos + d_k_matrix); 
            [ ~, fval ] = cplexlp(f_oos,Ain_oos,bin_oos); 
            fval_oos_dro_exantereg(i_set,i_eps) = fval;    
        end
    end
    oos_dro_exantereg(i_simu,:) = mean(fval_oos_dro_exantereg);
end
oos_reg_dro_exantereg = mean(oos_dro_exantereg) - oos_true;
    
% dro_expostreg
fval_oos_dro_expostreg = zeros(N_true_set,n_eps);
oos_dro_expostreg = zeros(n_simu,n_eps);
for i_simu = 1 : n_simu
    for i_set = 1 : N_true_set
        xi_hat_matrix_true = reshape( xi_vector_true_reshape(:,i_set), m, N_true_base )';
        dum_cell_x_true = cell(N_true_base, 1);
        for i_base = 1 : N_true_base
            dum_cell_x_true{i_base} = kron(eye(K), xi_hat_matrix_true(i_base,:)); 
        end    
        for i_eps = 1 : n_eps     
            x_oos = sol_dro_expostreg{i_simu,i_eps};        
            bin_oos = -blkdiag(dum_cell_x_true{:})*kron(ones(N_true_base,1),b_k_matrix) - kron(ones(N_true_base,1), c_k_trans_matrix*x_oos + d_k_matrix); 
            [ ~, fval ] = cplexlp(f_oos,Ain_oos,bin_oos); 
            fval_oos_dro_expostreg(i_set,i_eps) = fval;       
        end
    end
    oos_dro_expostreg(i_simu,:) = mean(fval_oos_dro_expostreg);        
end
oos_reg_dro_expostreg = mean(oos_dro_expostreg) - oos_true;

% dro_cost
fval_oos_dro_cost = zeros(N_true_set,n_eps);
oos_dro_cost = zeros(n_simu,n_eps);
for i_simu = 1 : n_simu         
    for i_set = 1 : N_true_set
        xi_hat_matrix_true = reshape( xi_vector_true_reshape(:,i_set), m, N_true_base )';
        dum_cell_x_true = cell(N_true_base, 1);
        for i_base = 1 : N_true_base
            dum_cell_x_true{i_base} = kron(eye(K), xi_hat_matrix_true(i_base,:)); 
        end
        for i_eps = 1 : n_eps 
            x_oos = sol_dro_cost{i_simu,i_eps};
            bin_oos = -blkdiag(dum_cell_x_true{:})*kron(ones(N_true_base,1),b_k_matrix) - kron(ones(N_true_base,1), c_k_trans_matrix*x_oos + d_k_matrix); 
            [ ~, fval ] = cplexlp(f_oos,Ain_oos,bin_oos); 
            fval_oos_dro_cost(i_set,i_eps) = fval; 
        end 
    end 
    oos_dro_cost(i_simu,:) = mean(fval_oos_dro_cost); 
end 
oos_reg_dro_cost = mean(oos_dro_cost) - oos_true;
