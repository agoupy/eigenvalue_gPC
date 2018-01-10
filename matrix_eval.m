function [M]=matrix_eval(points)
% Evaluate the matrix at the quadrature points

Nq=size(points,1);
M=cell(Nq,1);
C=[2 15 -8; 7 8 3; 233 0 -5];
A=[1 0 0; 0 0 0; 0 0 1];
B=[1 0 1; 0 1 0; 0 0 0];

for i=1:Nq
    M{i}=C+A.*points(i,1)+B.*points(i,2);
end

end