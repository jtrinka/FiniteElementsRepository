%% 2D Poisson Finite Element Approximation with Uniform Triangularation on a Square Domain

%Set up Domain
clear; clc;
xmin = 0;
xmax = 1;
ymin = xmin;
ymax = xmax;
numintpt = 20;
x = linspace(xmin,xmax,numintpt);
dx = x(2)-x(1);
y = linspace(ymin,ymax,numintpt);
dy = dx;
[x,y] = meshgrid(x,y);
x=x';
y=y';
% Define f(x,y)
f = @(x,y) 20*exp(-((x-.5).^2+(y-.5).^2)/0.05);

%Build the coefficient matrix
maindiag = 4*ones(numintpt-2,1);
A = diag(maindiag);
for j=1:numintpt-3
   A(j,j+1) = -1;
   A(j+1,j) = -1;
end

% I uncommented this block.  You were missing the banded structure of the
% matrix, but I'm not sure that this is the proper setup entirely.
% In fact, I'm sure this is incorrect.
for j=1:numintpt-5
   A(j,j+3) = -1;
   A(j+3,j) = -1;
end

%create basis functions
t = linspace(xmin,xmax,numintpt);
w = linspace(ymin,ymax,numintpt);
counter1=1; % I ADDED THIS IN
counter2=1;
for i=2:1:numintpt-1
    for j=2:1:numintpt-1
        
        %Region 1
        PQ1 = [x(i-1,j)-x(i,j);y(i-1,j)-y(i,j);-1];
        PR1 = [x(i,j-1)-x(i,j);y(i,j-1)-y(i,j);-1];
        Cross1 = cross(PQ1,PR1);
        R1 = @(t,w) ((((-Cross1(1,1)*(t-x(i,j))-Cross1(2,1)*(w-y(i,j)))/Cross1(3,1)))+1)*f(x(i,j),y(i,j));
        
        %Region 2
        PQ2 = [x(i-1,j)-x(i,j);y(i-1,j)-y(i,j);-1];
        PR2 = [x(i-1,j+1)-x(i,j);y(i-1,j+1)-y(i,j);-1];
        Cross2 = cross(PQ2,PR2);
        R2 = @(t,w) ((((-Cross2(1,1)*(t-x(i,j))-Cross2(2,1)*(w-y(i,j)))/Cross2(3,1)))+1)*f(x(i,j),y(i,j));
        
        %Region 3
        PQ3 = [x(i-1,j+1)-x(i,j);y(i-1,j+1)-y(i,j);-1];
        PR3 = [x(i,j+1)-x(i,j);y(i,j+1)-y(i,j);-1];
        Cross3 = cross(PQ3,PR3);
        R3 = @(t,w) ((((-Cross3(1,1)*(t-x(i,j))-Cross3(2,1)*(w-y(i,j)))/Cross3(3,1)))+1)*f(x(i,j),y(i,j));
        
        %Region 4
        PQ4 = [x(i,j+1)-x(i,j);y(i,j+1)-y(i,j);-1];
        PR4 = [x(i+1,j)-x(i,j);y(i+1,j)-y(i,j);-1];
        Cross4 = cross(PQ4,PR4);
        R4 = @(t,w) ((((-Cross4(1,1)*(t-x(i,j))-Cross4(2,1)*(w-y(i,j)))/Cross4(3,1)))+1)*f(x(i,j),y(i,j));
        
        %Region 5
        PQ5 = [x(i+1,j)-x(i,j);y(i+1,j)-y(i,j);-1];
        PR5 = [x(i+1,j-1)-x(i,j);y(i+1,j-1)-y(i,j);-1];
        Cross5 = cross(PQ5,PR5);
        R5 = @(t,w) ((((-Cross5(1,1)*(t-x(i,j))-Cross5(2,1)*(w-y(i,j)))/Cross5(3,1)))+1)*f(x(i,j),y(i,j));
        
        %Region 6
        PQ6 = [x(i+1,j-1)-x(i,j);y(i+1,j-1)-y(i,j);-1];
        PR6 = [x(i,j-1)-x(i,j);y(i,j-1)-y(i,j);-1];
        Cross6 = cross(PQ6,PR6);
        R6 = @(t,w) ((((-Cross6(1,1)*(t-x(i,j))-Cross6(2,1)*(w-y(i,j)))/Cross6(3,1)))+1)*f(x(i,j),y(i,j));
        
        %Region 1 line
        m1 = (y(i,j-1)-y(i-1,j))/(x(i,j-1)-x(i-1,j));
        b1 = (y(i-1,j)-m1*(x(i-1,j)));
        R1line = @(t) m1*t+b1;
        
        %Region 2/3 line
        m23 = (y(i,j)-y(i-1,j+1))/(x(i,j)-x(i-1,j+1));
        b23 = (y(i-1,j+1)-m23*(x(i-1,j+1)));
        R23line = @(t) m23*t+b23;
        
        
        %Region 4 line
        m4 = (y(i+1,j)-y(i,j+1))/(x(i+1,j)-x(i,j+1));
        b4 = (y(i,j+1)-m4*(x(i,j+1)));
        R4line = @(t) m4*t+b4;
        
        %Region 5/6 line
        m56 = (y(i+1,j-1)-y(i,j))/(x(i+1,j-1)-x(i,j));
        b56 = (y(i,j)-m56*(x(i,j)));
        R56line = @(t) m56*t+b56;
        
        %Region 1 integral
        R1int = integral2(R1,x(i-1,j),x(i,j),R1line,y(i,j));
        
        %Region 2 integral
        R2int = integral2(R2,x(i-1,j),x(i,j),y(i,j),R23line);
        
        %Region 3 integral
        R3int = integral2(R3,x(i-1,j),x(i,j),R23line,y(i,j+1));
        
        %Region 4 integral
        R4int = integral2(R4,x(i,j),x(i+1,j),y(i,j),R4line);
        
        %Region 5 integral
        R5int = integral2(R5,x(i,j),x(i+1,j),R56line,y(i,j));
        
        %Region 6 integral
        R6int = integral2(R6,x(i,j),x(i+1,j),y(i,j-1),R56line);
      
        
        H=R1int+R2int+R3int+R4int+R5int+R6int;
        
        RHS(counter1,:)=H;
        counter1 = counter1+1;
        counter2 = counter2+1;
    end
end

%Solve for interior solution
RHS = reshape(RHS,numintpt-2,numintpt-2);
uint = A\RHS;
% ********** This is the bulk of the issue ************
% The matrix A is the wrong size and RHS should be a vector (not a matrix).
%  The MATLAB backslash command is not doing what you think it is doing
%  here.  For numintpts = 20 we should have a 400x400 matrix A and a 400x1
%  vector for RHS.  In other words, I think that you're thinking of this as
%  a 1D problem but remember that we need to parse out the nodes in a 2D
%  problem so that we can use a matrix solve.

%Boundaries-Dirichlet
bound1 = zeros(numintpt-2,1);
utot1 = [bound1,uint,bound1];
bound2 = zeros(1,numintpt);
utot2 = [bound2;utot1;bound2];

%plot it
subplot(1,2,1)
surf(x,y,utot2) % this is a plot of the solution
subplot(1,2,2)
surf(x,y,f(x,y)) % this is a plot of the right-hand side: div(grad(u)) = -f
