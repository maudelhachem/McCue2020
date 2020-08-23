% Filename: PhasePlane.m
% Author: Maud El-Hachem
% School of Mathematical Sciences, 
% Queensland University of Technology, Brisbane, Australia.

% Reference: Scott W McCue, Maud El-Hachem, Matthew J Simpson (2020)
% Exact sharp-fronted travelling wave solutions.

% This function solves numerically Equations (11) and (12) 
% in the phase plane by Heun's method and displays two travelling wave 
% solutions in the phase plane such as U(z)>=0.
% The same figure shows also the equilibrium points (0,0) and (1,0).
% The invading travelling wave is obtained from integrating Equations (11)
% and (12) in the positive direction of z, starting on unstable manifold of
% the saddle (1,0).
% The receding travelling wave is obtained from integrating Equations (11)
% and (12) in the negative direction of z, starting on stable manifold of
% the saddle (1,0).

% Wave speed c
c = 5/sqrt(6);

% Calculating the eigenvalues and the eigenvectors of the solution 
% around the saddle (1,0)
Us = 1;
A = [0 1;(2*Us-1) -c];
[v,d] = eig(A);

% Determining the initial conditions close to the equilibrium point
% (1,0) along the eigenvector of the solution
% Setting the initial conditions to start on the unstable manifold
% To obtain the invading travelling wave solution
ICinv = [0;0];
if (d(1,1) > 0)
    ICinv = v(:,1);
end
if (d(2,2) > 0)
    ICinv = v(:,2);
end
% Initial conditions to obtain the invading travelling wave solution
ICinv = [Us+ICinv(2,1)*0.001; ICinv(1,1)*0.001;];

% Setting the initial conditions to start on the stable manifold
ICrec = [0;0];
if (d(1,1) < 0)
    ICrec = v(:,1);
end
if (d(2,2) < 0)
    ICrec = v(:,2);
end
% Initial conditions to obtain the receding travelling wave solution
ICrec = [Us+ICrec(1,1)*0.001; ICrec(2,1)*0.001;];

% Solving Equations Equations (11) and (12) in the phase plane with
% Heun's method
% Step size dz used to discretise the domain of z
dz = 0.001;
% Lower and upper limit of the domain of z such as z_begin <= z <= z_end
% The initial conditions are applied at z = z_begin.
z_begin = 0;
z_end = 30;
% Solving for the invading travelling wave by integrating in the positive
% direction of z from z_begin to z_end
% The solution range is limited to [0 3] for U and [-1 3] for V.
[Uinv, Vinv] = heunSolver(c, dz, z_begin, z_end, ICinv, [0 2], [-1 3]);
% Solving for the receding travelling wave by integrating in the negative
% direction of z from z_begin to -z_end
% The solution range is limited to [0 3] for U and [-1 3] for V.
[Urec, Vrec] = heunSolver(c, -dz, z_begin, -z_end, ICrec, [0 3], [-1 3]);

% Setting and computing the field vectors of the solution
xmin = -1;
xmax = 3;
ymin = -1.0;
ymax = 3;
y1 = linspace(xmin,xmax,14);
y2 = linspace(ymin,ymax,15);
[x,y] = meshgrid(y1,y2);
du = zeros(size(x));
dv = zeros(size(x));
f1 = @(Y) [Y(2);(-c*Y(2)+(-Y(1)+Y(1).*Y(1)));];
for i = 1:numel(x)
    Yprime = f1([x(i); y(i)]);
    du(i) = Yprime(1);
    dv(i) = Yprime(2);
end

% Displaying black lines to indicate the axis at U(z)=0 and V(z)=0
line([xmin xmax],[0 0],'Color','k','LineStyle','-','LineWidth',0.5);
hold on
line([0 0],[ymin ymax],'Color','k','LineStyle','-','LineWidth',0.5);
% Displaying the travelling waves solutions  
plot(Uinv,Vinv,'r-','LineWidth',3);
plot(Urec,Vrec,'r-','LineWidth',3);
% Displaying the field vectors
q = quiver(x,y,du*0.05,dv*0.05,'b','LineWidth',1,'AutoScale','on'); 
q.ShowArrowHead = 'off';
% Displaying the equilibrium points (0,0) and (1,0)
plot(0,0,'o','MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',4);
plot(1,0,'o','MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',4);
% Displaying the labels of the corresponding axis
ylabel('$V(z)$','interpreter','latex');
xlabel('$U(z)$','interpreter','latex'); 
% Setting the limits of the current axis U(z) and V(z)
xlim([xmin,xmax]);
ylim([ymin,ymax]);

% Setting fonts
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(gca,'fontsize', 18);
box on;
hold off


% Function heunSolver
% This function solves Equations (10) and (11) by Heun's method 
% INPUT ARGUMENTS:
% ** c, the wave speed, if positive, 0 < c < 2, if negative, c < 0 
% ** dz, the step size used to discretise the domain of z
% ** z_begin and z_end, the lower and upper limit of the numerical domain 
% of z such as z_begin <= z <= z_end. The initial conditions are applied at
% z = z_begin.
% ** IC = [U0, V0], the values of the initial conditions
% ** Ulimits: [Umin Umax] the boundaries of the range for U(z)
% ** Vlimits: [Vmin Vmax] the boundaries of the range for V(z)
% OUTPUT ARGUMENTS:
% ** Uout : The solution U(z)
% ** Vout : The solution V(z)
% The output size of Vout and Xout may not be equal to the original number
% of nodes correponding to (z_end-z_begin)/dz+1.
% The array Uout and Vout may be truncated. 

function [Uout, Vout] = heunSolver(c,dz,z_begin,z_end,IC,Ulimits,Vlimits)

    % Domain of z
    z = z_begin:dz:z_end;
    % Number of nodes in the domain
    sz = length(z);

    % Initialisation 
    V = zeros(sz,1);
    U = zeros(sz,1);
    U(1) = IC(1);
    V(1) = IC(2);
   
    % Current output size of the domain
    szout = sz;

    % Integrating Equations (11) and (12) with Heun's method
    for i = 1:sz-1
        Ubar = dz * V(i) + U(i); 
        Vbar = dz * (-c*V(i)- U(i)*(1-U(i))) + V(i);
        U(i+1) = dz/2 * (V(i)+Vbar) + U(i); 
        V(i+1) = dz/2 * ((-c*V(i) - U(i)*(1-U(i))) + (-c*Vbar - Ubar *(1-Ubar))) + V(i); 
        % If the solutions are outside the required ranges for U and V, stop the solver
        if (V(i+1) < Vlimits(1) || V(i+1) > Vlimits(2) || U(i+1) < Ulimits(1) || U(i+1) > Ulimits(2))
            szout = i;
            break;
        end
    end

    % Truncating the size of the output vectors
    Uout = U(1:szout,1);
    Vout = V(1:szout,1);
end