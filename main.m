%% gPC decomposition of the eigenvalues and eigenvectors of a random matrix

% Based on 'Efficient characterization of the random eigenvalue problem 
% in a polynomial chaos decomposition' by Roger Ghanem and Debraj Ghosh

%% gPC decomposition of the matrix
% The random matrix has to be defined in matrix_eval

order_gPC_matrix=3;                                                         % gPC order for matrix decomposition
dim=2;                                                                      % number of random parameter
K=matrix_gPC(order_gPC_matrix,dim);
L=size(K,1);

%% Construction of relevant matrices B and Gamma
order_gPC_eigen=3;
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

return;
%% Evaluate gPC convergence with quasi Monte-Carlo comparison
N_QMC=500;
points=Sobol_quasi_random(N_QMC,dim);
ev_QMC=zeros(3,N_QMC);
M=matrix_eval(points);

for i_QMC=1:N_QMC
    ev_QMC(:,i_QMC)=eig(M{i_QMC})';
end