function [ y ] = f_test( x )
%function for integration testing

%y=(-x+1).*(x>=0).*(x<=1)+(x+1).*(x<0).*(x>=-1);

%y=cos(pi*x/2).*(x<=1).*(x>=-1);

y=exp(-3*x.^2);

%y=exp(-3*x(:,1).^2)+exp(-(x(:,2)-1).^2);

end

