function [M] = assembly(A,b)
m = size(A,1)+size(b,1)-1;
n = size(A,2)+size(b,2)-1;
M1 = zeros(m,n);
M2 = zeros(m,n);
M1(1:size(A,1),1:size(A,2))=A;
M2(size(A,1):m,size(A,2):n)=b;
M = M1+M2; 
end
