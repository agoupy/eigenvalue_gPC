%% gPC decomposition of the eigenvalues and eigenvectors of a random matrix

% Based on 'Efficient characterization of the random eigenvalue problem 
% in a polynomial chaos decomposition' by Roger Ghanem and Debraj Ghosh

%% gPC decomposition of the matrix
% The random matrix has to be defined in matrix_eval

order_gPC_matrix=4;                                                         % gPC order for matrix decomposition
dim=2;                                                                      % number of random parameter
K=matrix_gPC(order_gPC_matrix,dim);
L=size(K,1);

%% Construction of relevant matrices B and Gamma
order_gPC_eigen=6;
alpha=multi_index(dim,order_gPC_eigen);
P=size(alpha,1);                                                            % number of polynomials

n=size(K{1},1);

B=cell(P,1);
Gamma=cell(P,1);

for i=1:L
    B{i}=kron(eye(P),K{i});
    Gamma{i}=construct_gamma(dim,order_gPC_eigen,n,i);
end

% Complete with null matrices /!\ P has to be greater than L
for i=L+1:P
    B{i}=zeros(n*P);
    Gamma{i}=construct_gamma(dim,order_gPC_eigen,n,i);
end

%% Find coefficients with Newton-Raphson approach

[ev,ef]=newton_raphson(B,Gamma,K{1});


%% Evaluate gPC convergence with quasi Monte-Carlo comparison
N_QMC=10000;
points=Sobol_quasi_random(N_QMC,dim);
ev_QMC=zeros(n,N_QMC);
M=matrix_eval(points);

for i_QMC=1:N_QMC
    ev_QMC(:,i_QMC)=eig(M{i_QMC})';
end

% Mean
ev_mean_QMC=mean(ev_QMC,2);
ev_mean_gPC=zeros(n,1);
for i=1:n
    ev_mean_gPC(i,1)=ev{i}(1,1);
end

% Variance
ev_var_QMC=var(ev_QMC,0,2);
ev_var_gPC=zeros(n,1);
for i=1:n
    ev_var_gPC(i,1)=sum(abs(ev{i}(1,2:end).^2));
end

display_result(ev_mean_gPC,ev_mean_QMC,ev_var_QMC,ev_var_gPC);