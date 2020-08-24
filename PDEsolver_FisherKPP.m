% Filename: PDEsolver_FisherKPP.m
% Author: Maud El-Hachem
% School of Mathematical Sciences, 
% Queensland University of Technology, Brisbane, Australia.

% Reference: Scott W McCue, Maud El-Hachem, Matthew J Simpson (2020)
% Exact sharp-fronted travelling wave solutions.

% This script solves the Fisher-KPP equation (1) with the IC conditions (2)
% and boundary conditions (3).
% It displays the figure showing the density profiles at time
% t = 0, 3, 6 and 9, as in Figure 2(a).

% Tolerance 
tol = 1e-08;
% Total time
total_time = 9;
% Time step dt
dt = 0.01;
% Total number of time steps
ts = round(total_time/dt+1);
% Times when to print the density profiles
time_toPrint = [0 3 6 9];
% Colors array used to print density profiles at required time steps
colors = [1.00 0.48 0; 0.0471 0.5098 0;0 1.00 1.00; 1.00 0 1.00; ];
% Spatial grid step dx
dx = 0.0001;
% Spatial domain
L = 60;
x = 0:dx:L;
% Number of nodes in the spatial domain
nodes = size(x,2);
% Initialisation of variables used in Newton-Raphson algorithm
% Correction of densities at each iteration
delta_u = ones(1,nodes);
% Function F
Fu = zeros(1,nodes);
% Densities at time j+1
u = zeros(1,nodes);

% Coefficients a b c of the tridiagonal matrix
% Jacobian
% J(u)
coeffA_u = zeros(1,nodes);
coeffB_u = zeros(1,nodes);
coeffC_u = zeros(1,nodes);

% Initialisation of density u(x,0)
for i = 1:nodes
    if (x(i) < 10)
        u(1,i) = 0.5;
    else
        u(1,i)=0.5*exp(-2/sqrt(6)*(x(i)-10));
    end
        
end
% Initialistion of the densities at time j
u_p = u;

% Opening a new figure
figh = figure;

hold on
% Displaying the initial conditions
if (time_toPrint(1)==0)
    plot(x(1:nodes), u, '-', 'LineWidth',2 ,'Color', colors(1,1:3));
end

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(gca,'fontsize', 18);
%% Newton-Raphson algorithm - main loop
% A variable to calculate the position where U = 0.5
xhalf = 0;
% For all time steps
for j = 1:ts
    % Calulating and displaying current time
    t = j * dt;
    disp(t);

    condition = 1;

    % While the correction term is not small enough
    while (condition)

        % Boundary condition at x = 0 as in Equation (3)
        % Supplementary Material Equation (S4)
        coeffA_u(1,1) = 0.0;
        coeffB_u(1,1) = -1.0;
        coeffC_u(1,1) = 1.0;
        Fu(1) = -1.0*(u(1,2)-u(1,1));

        % Boundary condition at x = L as in Equation (3)
        % Supplementary Material Equation (S5)
        coeffA_u(1,nodes) = 0.0;
        coeffB_u(1,nodes) = 1.0;
        coeffC_u(1,nodes) = 0.0;
        Fu(1,nodes) = -1.0*(u(1,nodes));

        % J(u) delta u = -F(u)
        for i = 2:nodes-1
           coeffA_u(1,i) = 1.0/(dx^2);
           coeffB_u(1,i) = - 2.0/(dx^2)  - 1.0/dt + (1 - 2.0*u(1,i));
           coeffC_u(1,i) = 1.0/(dx^2);
           Fu(1,i) = -((u(1,i+1) - 2*u(1,i) + u(1,i-1))/(dx^2) ...
               - (u(1,i)-u_p(1,i)) / dt + u(1,i)*(1-u(1,i)));
        end 
        delta_u = tridia(coeffA_u, coeffB_u, coeffC_u, Fu, nodes);

        % Correction of u
        for i = 1:nodes
            u(1,i) = u(1,i) + delta_u(1,i);
        end
        % Checking if \Delta u is small enough
        % \Delta u <= tolerance
        if (norm(delta_u,Inf) <= tol)
            condition = 0;
        end

    end

    % Displaying density profiles solutions at required times
    if (isempty(find(time_toPrint == t,1)) == 0)
        textc = strcat(strcat('$t = ',num2str(t)), '$');
        % Choosing a different color for the last required 
        if (t == time_toPrint(end))
            colorNo = 3;
        else
            colorNo = 2;
        end
        plot(x(1:nodes), u, '-', ...
            'LineWidth',2, 'DisplayName', textc, 'color', ...
            colors(colorNo,1:3));
    end
    
    % Calculating the wave speed as described in the Supplementary Material
    % section S2
    for i = 1:nodes
        if (u(1,i) > 0.5 && u(1,i+1) < 0.5)
            c = (x(i) - xhalf)/dt;
            xhalf = x(i);
            break;
        end
    end

    % Updating current u(x,t) 
    u_p = u;
end

%%
% Displaying the exact solution shifted such as U(z) = 0.5 at z=0
xexact=-40:dx:40;
solutionExact = 1./(1+(sqrt(2)-1)*exp(xexact/sqrt(6))).^2;
plot(xexact+xhalf,solutionExact,':','LineWidth',2,'Color', colors(4,1:3));

% Setting the limit of the current axis
xlim ([0 40])

% Displaying axis labels
ylabel('$u(x,t)$','interpreter','latex','fontsize',18);
xlabel('$x$','interpreter','latex','fontsize',18);

% Setting fonts
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(gca,'fontsize', 18);
box on

hold off
%% Function tridia
% This function implements Thomas algorithm that solves a system
% of equations Ax = d, where A is a tridiagonal matrix. The parameters 
% a,b and c are the three diagonals of the matrix A. N is the size of 
% the vector solution x.
function x = tridia(a,b,c,d,N)
    x=zeros(1,N);
    bb=b;
    dd=d;
    for i=2:N
        ff=a(i)/bb(i-1);
        bb(i)=bb(i)-c(i-1)*ff;
        dd(i)=dd(i)-dd(i-1)*ff;
    end

    for i=1:N-1
    x(N)=dd(N)/bb(N);    
    j=N-i;
    x(j)=(dd(j)-c(j)*x(j+1))/bb(j);
    end
end 
