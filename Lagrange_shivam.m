function [ Sig ] = Lagrange_shivam(n,start,he)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes heren=

%This function creates the lagrange interpolation in general for any order
%and for any start and end point 


%n=input('enter the order:');
%start=input('enter Start Point:');
%he=input('enter end point:');


P=1;
syms x

for i=start+1:n+1
    for j=start+1:n+1
        if (j-1)~=(i-1)
        P=P*(x-(j-1)*he/n)/((i-1)*he/n -(j-1)*he/n);
        end
    end
    Sig(i)=P;
    P=1;
end

 
end




