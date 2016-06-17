%% One-Spatial Dimension Poisson Finite Elements
% Sullivan: To check this code I've used a right-hand side that has an
% easy-to-find analytic solution: f(x) = x(1-x).  The solution is 
% u(x) = -x^3/6 + x^4/12 + (1/6-1/12)x.  You can check that if you take the
% negative second derivative of this u(x) you will get f(x).  I added a
% plot to show that your solution gives an approximation to the analytic
% solution.  I also added an error plot since I have the analytic solution.
%
% There were some mistakes 
clear; clc;
xmin=0;
xmax=1;
npt=50; 
% WITH A SMALL NUMBER OF GRID POINTS I FOUND THAT YOU WERE INCORRECTLY 
% APPROXIMATING THE LEFT-HAND SIDE.  WITH A LARGE NUMBER OF GRID POINTS
% THIS MISTAKE WAS INVISIBLE.
% YOU CAN ALSO SEE THAT EVEN WITH AS FEW AS 50 INTERIOR POINTS THE MAXIMUM
% ERROR IS ON THE ORDER OF 10^(-6).  TAKING 10,000 POINTS IS TOTAL
% OVERKILL.

x=linspace(xmin,xmax,npt)';
f = @(x) x*(1-x); % this is a simple function where I can check the answer
g = @(x) -x.^3/6 + x.^4/12 +(1/6-1/12)*x; % this is the analytic solution
% f = @(x) 5*sin(2*pi*x.^10); % OLD
boundary=0;
%uniform partitions
dx= x(2)-x(1);
%create basis functions
t = linspace(xmin,xmax,npt)';



counter=1; % I ADDED THIS IN
for j = 2:1:length(x)-1 % YOU HAD A MISTAKE HERE
    prodfplus  = @(t) ((t-x(j-1))/dx)*f(x(j));
    prodfminus = @(t) ((-t+x(j+1))/dx)*f(x(j));
    
    F=SimpsonRulef(prodfplus,x(j-1),x(j),length(x));
    
    G=SimpsonRulef(prodfminus,x(j),x(j+1),length(x));
     
    H= (F+G);
      
    b(counter,:)=H; % I'M INCREMENTING BY THE COUNTER HERE
    counter=counter+1; % INCREMENT THE COUNTER
% I FIXED THE CRAZY INDENTATIONS THAT YOU HAD ... CLEAN CODE = :-)      
    
end




maindiag = 2*ones(npt-2,1);
A = diag(maindiag);
for j=1:npt-3
   A(j,j+1) = -1;
   A(j+1,j) = -1;
end

A = (1/dx)*A;
uin= A\b; %coefficeint backslash vector
u = [0;uin;0];
subplot(2,1,1)
plot(x,u,'b',x,g(x),'r--')
subplot(2,1,2)
plot(x,abs(g(x)-u),'k--')