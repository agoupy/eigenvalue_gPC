%% gPC decomposition of the eigenvalues and eigenvectors of a random matrix

% Based on 'Efficient characterization of the random eigenvalue problem 
% in a polynomial chaos decomposition' by Roger Ghanem and Debraj Ghosh

%% gPC decomposition of the matrix
% The random matrix has to be defined in matrix_eval

L=4;                        % gPC order for matrix decomposition
dim=2;                      % number of random parameter
K=matrix_gPC(L,dim);

%% Construction of relevant matrices B and Gamma


%% Find coefficients with Newton-Raphson approach