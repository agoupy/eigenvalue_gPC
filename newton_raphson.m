function [ev,ef]=newton_raphson(B,Gamma,K)
% Compute the gPC coefficients by a newton raphson procedure

P=size(B,1);
n=size(K,1);
ev=cell(n,1);
ef=cell(n,1);
% collect expectations in Gamma
expectation=zeros(P,P,P);
for i=1:P
    for j=i:P
        for k=j:P
            expect=Gamma{i}(n*(j-1)+1,n*(k-1)+1);
            expectation(i,j,k)=expect;
            expectation(i,k,j)=expect;
            expectation(j,i,k)=expect;
            expectation(j,k,i)=expect;
            expectation(k,j,i)=expect;
            expectation(k,i,j)=expect;
        end
    end
end


%% Initial points for the optimization
[ef_init,ev_init]=eig(K);
for i_ev=1:n
    aux_ev=zeros(1,P);
    aux_ev(1,1)=ev_init(i_ev,i_ev);
    ev{i_ev}=aux_ev;
    aux_ef=zeros(n,P);
    aux_ef(:,1)=ef_init(:,i_ev);
    ef{i_ev}=aux_ef;
end

%% NR iterations

for i_ev=1:n
    X=zeros((n+1)*P,1);
    X(1,1)=ev{i_ev}(1,1);
    X(P+1:P+n,1)=ef{i_ev}(:,1);
    
    F=func_eval(X,B,Gamma);
    J=jacob(X,B,expectation);
    tol=1e-6;
    iter=0;
    
    while (max(abs(F))>tol) && (iter<1000)
        X=X-J\F;
        F=func_eval(X,B,Gamma);
        J=jacob(X,B,expectation);
        iter=iter+1;
    end
    
    if iter==1000
        disp('Max iteration reached ! /!\ NR not converged')
    else
        disp(['Zero reached in ' num2str(iter) ' iterations'])
    end
    
    ev{i_ev}=X(1:P,1)';
    ef{i_ev}=reshape(X((P+1):end,1),n,P);
end
end


%% Function evaluation
function F=func_eval(X,B,Gamma)
% function evaluation F(X)

P=size(B,1);
n=size(X,1)/P-1;

F=zeros((n+1)*P,1);
AA=zeros(n*P);
phi=X((P+1):end,1);

for ii=1:P
    AA=AA+B{ii}*Gamma{ii}-X(ii,1)*Gamma{ii};
end

F(1:n*P,1)=AA*phi;

for kk=1:P
    F(n*P+kk,1)=phi'*Gamma{kk}*phi-1*(kk==1);
end

end


%% Computation of the jacobian
function J=jacob(X,B,expectation)
% compute the jacobian of the former function

P=size(B,1);
n=size(X,1)/P-1;
phi=X((P+1):end,1);
J=zeros((n+1)*P,(n+1)*P);

% eq 37
for ii=1:P
    for kk=1:P
        sum=zeros(n,1);
        for jj=1:P
            sum=sum-expectation(ii,jj,kk)*phi(((jj-1)*n+1):(jj*n));
        end
        J(((kk-1)*n+1):(kk*n),ii)=sum;
    end    
end

% eq 38
for kk=1:P
    for jj=1:P
        sum=zeros(n);
        for ii=1:P
            sum=sum+expectation(ii,jj,kk)*(B{ii}(1:n,1:n)-X(ii,1)*eye(n));
        end
        J(((kk-1)*n+1):(kk*n),(P+(jj-1)*n+1):(P+jj*n))=sum;
    end
end

% eq 39
% Nothing to do :) !

% eq 40
for kk=1:P
    for jj=1:P
        for ll=1:n
            sum=0;
            for ii=1:P
                sum=sum+2*phi((ii-1)*n+ll)*expectation(ii,jj,kk);
            end
            J(n*P+kk,P+(jj-1)*n+ll)=sum;
        end
    end
end


end