function [K]=matrix_gPC(L,dim)
% Compute the gPC decomposition of the random matrix using quadrature
% The matrix is evaluated in matrix_eval


%% Create quadrature
addpath('/InterfaceMATLAB/');
domain = repmat([ 0, 1/2],[dim,1]);                                         % used to have exp(-x^2/2) for prob and exp(-x^2) for phys
[ weights, points ] = tsgMakeQuadrature( dim,'gauss-hermite','qptotal',L,0,domain,0);
weights=weights./((2*pi)^(dim/2));                                          % normalization pi for phys 2*pi for prob
Nq=size(points,1);

%% Call matrix_eval
M=matrix_eval(points);

%% Compute gPC coefficients
%generate the hermite polynomials
alpha=multi_index(dim,L);
He=poly1D(L,'hermite-prob');                                                % Hermite polynomials of degree 1
P=size(alpha,1);                                                            % number of polynomials

%compute the coefficients
K=cell(P,1);
for i=1:P
    %polynomial evaluation
    pol=ones(Nq,1);
    for j=1:dim
        pol=polyval(He{alpha(i,j)+1},points(:,j)).*pol;
    end
    %normalisation for phys
    %pol=pol./sqrt(prod(factorial(alpha(num,:)).*2.^alpha(num,:)));
    %normalisation for prob
    pol=pol/sqrt(prod(factorial(alpha(i,:))));
    
    %integration
    K{i}=0;
    for kk=1:Nq
        K{i,1}=K{i,1}+M{kk,1}*pol(kk,1)*weights(kk,1);
    end
end


end