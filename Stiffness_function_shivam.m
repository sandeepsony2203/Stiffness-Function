
function [K,Assembled_K,f] = Stiffness_function_shivam()
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

%Example :If we devide the system i.e length 0 to L into 4 elements and
%each single element is lineary divided.

% enter the number elements in which the system is being divided :4
% Division b/w any two elements1
% enter Start Point:0
% enter end point of your system:1
% enter the order2

%n=input('enter the dimension for which stiffness matrix needs to be calculated:');
clear all;
clc;
e=input('enter the number elements in which the system is being divided :');
n=input('Division b/w any two elements :');
start=input('enter Start Point :');
last=input('enter end point of your system :');
order=n+1;%input('enter the order');

a=1;
c=1;

syms x y 
L=Lagrange_shivam(n,start,last/e);
%F=@(x,k)(-(x+((k-1)*last/n))^2); % change here for F
F=@(x,k)(-(x+((k-1)*last/e))^2);
Assembled_K = [];
for i=1:n+1
    for j=1:n+1
        
        %P(i,j)=(int(a*diff(L(1,i))*diff(L(1,j))+c*L(1,i)*L(1,j) , start ,last)); % change here for K
        P(i,j)=(int(diff(L(1,i))*diff(L(1,j))-L(1,i)*L(1,j) , start ,last/e));
        Integ=a*diff(L(1,i))*diff(L(1,j))-c*L(1,i)*L(1,j);
        K(i,j)=Gauss_Quadrature_interation_shivam(Integ,start,last/e,order);
        
    end
    
end

for i=1:e
    for j=1:n+1
        f(i,j)=int(F(x,i)*L(1,j),start,last/e);
    end
end


 Assembled_K = P;
 for i=2:e
  Assembled_K = assembly(Assembled_K,P);
 end
 
 disp('The Assembled matrix K :')
 Assembled_K
  Assembled_F = (f(1,:))';
 for i=2:e
  Assembled_F = assembly(Assembled_F,(f(i,:))');
 end
 disp('The Assembled matrix F :')
 Assembled_F
 
 u=zeros(n+1,1);
 u(1)=input('enter initial and final boundary[ u_i ]  :');
 u(e+1)=input('enter initial and final boundary[ u_i] :');
 
 Q=input('enter Q :');
 
 y=  inv(Assembled_K(2:e,2:e))*(Assembled_F(2:e,1)+Q(2:e,1));
 u(2:e,1)=y(1:end,1);
 disp('Finally:')
 disp('Final U:')
 u
 disp('Final Q :')
 D=Assembled_K*u-Assembled_F
 
 plot(u,'r')
 title('Calculted U')
 hold on 
 
%  gt=@(l)(-0.91*sin(l) +0.5*cos(l)+l^2-1/2);
%  for i=1:e+1
%      v(i)=gt(start +(e-1)*last/e);
%  end
% plot(v,0,'rs')
%      hold on
end
