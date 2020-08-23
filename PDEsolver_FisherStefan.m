% Filename: FisherStefan.m
% Author: Maud El-Hachem

% School of Mathematical Sciences, 
% Queensland University of Technology, Brisbane, Australia.

% Reference: Scott W McCue, Maud El-Hachem, Matthew J Simpson (2020)
% Exact sharp-fronted travelling wave solutions.

% This script solves the Fisher-Stefan equation (4) 
% with boundary conditions (5)-(6).
% It displays the figure showing the density profiles at time
% t = 0, 3, 6 and 9, as in Figure 2(b).

% Parameter kappa of the Stefan condition that is associated to a wave
% speed c = -5/sqrt(6)
kappa = -0.906610965581149;
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
% Spatial grid step dxi
dxi = 0.00001;
% Initial position of the moving boundary
st = 100;
% Initialisation the moving boundary at time step j
st_p = st;
st_array = zeros(1,ts); 
% Spatial domain as in Supplementary Material Section S4
xi = 0:dxi:1;
% Number of nodes in the spatial domain
nodes_xi = size(xi,2);

% Initialisation of variables used in Newton-Raphson algorithm
% Correction of densities at each iteration
delta_u = ones(1,nodes_xi);
% Function F
Fu = zeros(1,nodes_xi);
% Initialistion of the densities at time j
u = zeros(1,nodes_xi);

% Coefficients a b c of the tridiagonal matrix
% Jacobian
% J(u)
coeffA_u = zeros(1,nodes_xi);
coeffB_u = zeros(1,nodes_xi);
coeffC_u = zeros(1,nodes_xi);

% Initialisation of density u(xi,0)
for i = 1:nodes_xi
    u(1,i) = 0.5;
end
u(1,nodes_xi) = 0;

% Initialistion of the densities at time j
u_p = u;

% Opening a new figure
figh = figure;

hold on
% Displaying the initial conditions
if (time_toPrint(1)==0)
    plot(xi(1:nodes_xi)*st, u, '-', 'LineWidth',2 ,...
        'DisplayName', '$u(x,0)$','Color', colors(1,1:3));
end

%% Newton-Raphson algorithm - main loop
% For all time steps
for j = 1:ts
    % Current time
    t = j * dt;
    disp(t);
    % Current position of the moving boundary
    st_array(j) = st;
    condition = 1;

    % While the correction term is not small enough
    while (condition)

        % Boundary condition at x = 0 from Equation (5)
        % Supplementary Material Equation (S14)
        coeffA_u(1,1) = 0.0;
        coeffB_u(1,1) = -1.0;
        coeffC_u(1,1) = 1.0;
        Fu(1) = -1.0*(u(1,2)-u(1,1));

        % Boundary condition at x = L from Equation (6)
        % Supplementary Material Equation (S15)
        coeffA_u(1,nodes_xi) = 0;
        coeffB_u(1,nodes_xi) = 1.0;
        coeffC_u(1,nodes_xi) = 0;
        Fu(1,nodes_xi) = -1.0 * (u(1,nodes_xi));

        % J(u) delta u = -F(u)
        for i = 2:nodes_xi-1
           coeffA_u(1,i) = 1.0/(dxi^2*st^2) - (i-1)*dxi/st * (st-st_p)/(2*dt*dxi);
           coeffB_u(1,i) = - 2.0/(dxi^2*st^2)  - 1.0/dt + (1 - 2.0*u(1,i));
           coeffC_u(1,i) = 1.0/(dxi^2*st^2) + (i-1)*dxi/st * (st-st_p)/(2*dt*dxi);
           Fu(1,i) = -(u(1,i+1) - 2*u(1,i) + u(1,i-1))/(dxi^2*st^2) ...
               - (i-1)*dxi/st * (u(1,i+1) - u(1,i-1)) * (st-st_p)/(2*dt*dxi) ...
               + (u(1,i)-u_p(1,i)) / dt - u(1,i)*(1-u(1,i));
        end 
        delta_u = tridia(coeffA_u, coeffB_u, coeffC_u, Fu, nodes_xi);

        % Correction of u
        for i = 1:nodes_xi
            u(1,i) = u(1,i) + delta_u(1,i);
        end
        
        % Checking if \Delta u is small enough
        % \Delta u <= tolerance
        if (norm(delta_u,Inf) <= tol)
            condition = 0;
        end
        
        % Stefan condition from Equation (6)
        % Supplementary Material Eq (S16)
        st = st_p + dt*(-kappa*(u(1,nodes_xi-2)/2-2*u(1,nodes_xi-1))/(dxi*st));
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
        plot(xi(1:nodes_xi)*st, u, ...
            'LineWidth',2, 'DisplayName', textc, 'color', colors(colorNo,1:3));
    end

    % Updating current s(t)
    st_p = st;
    % Updating current u(x,t)
    u_p = u;

end
    
% Calculating the wave speed as described in the Supplementary Material
% section S3
c = (st_array(end) - st_array(end-1))/dt;
for i = 1:nodes_xi
    if (u(1,i) > 0.5 && u(1,i+1) < 0.5)
        xhalf = xi(i);
        break;
    end
end
%%

% Setting the limit of the current axis
xlim([60 100])
ylim([0 1])

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
