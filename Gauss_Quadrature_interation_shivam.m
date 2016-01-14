function [ S ] = Gauss_Quadrature_interation_shivam(Integ,start,last,order)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here


syms t x

c1=(start+last)/2;
c2=(last-start)/2;

S=0;
if order==1
     w(1)=0;
     W(1)=2;
     s=matlabFunction(c2*subs(Integ,x,c1+c2*t));
elseif order==2
        w(1)= -1/sqrt(3);w(2)= 1/sqrt(3);
        W(1)=1;W(2)=1;
        subs(Integ,x,c1+c2*t);
        s=matlabFunction(c2*subs(Integ,x,c1+c2*t));
elseif order==3
        w(1)=5/9; w(2)=8/9; w(3)=5/9;
        W(1)=-sqrt(3/5);W(2)=0;W(3)=sqrt(3/5);
         s=matlabFunction(c2*subs(Integ,x,c1+c2*t));
end

        for i=1:order 
            S=S+W(i)*s(w(i));
         end

end
