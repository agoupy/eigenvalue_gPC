function [M]=matrix_eval(points)

dim=size(points,2);
Nq=size(points,1);
M=cell(Nq,1);
C=[2 3 4; 5 8 2; 7 5 3];
A=[0 1 0; 0 0 1; 1 0 0];
B=[1 0 0; 1 0 0; 0 0 1];

for i=1:Nq
    M{i}=A.*points(i,1)+B.*points(i,2)+C;
end

end