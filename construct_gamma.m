function Gamma=construct_gamma(dim,order,n,i)
% Build the Gamma matrix n°i using Hermite polynomials

%% Create quadrature
addpath('~/eigenvalue_gPC/InterfaceMATLAB/');
domain = repmat([ 0, 1/2],[dim,1]);                                         % used to have exp(-x^2/2) for prob and exp(-x^2) for phys
[ weights, points ] = tsgMakeQuadrature(dim,'gauss-hermite','qptotal',dim*order,0,domain,0);
weights=weights./((2*pi)^(dim/2));                                          % normalization pi for phys 2*pi for prob
Nq=size(points,1);


%% Polynomial definition
alpha=multi_index(dim,order);
P=size(alpha,1);                                                           % number of polynomials
He=poly1D(order,'hermite-prob');                                           % Hermite polynomials of degree 1
pol=cell(P,1);
Gamma=zeros(n*P);
ortho=zeros(P);

% Compute polynomials on quadrature points
for i_poly=1:P
    pol{i_poly,1}=ones(Nq,1);
    for j=1:dim
        pol{i_poly,1}=polyval(He{alpha(i_poly,j)+1},points(:,j)).*pol{i_poly,1};
    end
    %normalisation for phys
    %pol=pol./sqrt(prod(factorial(alpha(num,:)).*2.^alpha(num,:)));
    %normalisation for prob
    pol{i_poly,1}=pol{i_poly,1}/sqrt(prod(factorial(alpha(i_poly,:))));
end


%% Compute expectation and fill in matrix
for j=1:P
    for k=j:P
        % Compute E[He{i}He{j}He{k}]
        expectation=(pol{i,1}.*pol{j,1}.*pol{k,1})'*weights;
        ortho(j,k)=(pol{j,1}.*pol{k,1})'*weights;
        ortho(k,j)=(pol{j,1}.*pol{k,1})'*weights;
            
        % Fill in Gamma
        Gamma((j-1)*n+1:j*n,(k-1)*n+1:k*n)=expectation*eye(n);
        Gamma((k-1)*n+1:k*n,(j-1)*n+1:j*n)=expectation*eye(n);
    end
end

end