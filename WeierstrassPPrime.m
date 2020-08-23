% Filename: WeierstrassPPrime.m
% Author: Maud El-Hachem
% School of Mathematical Sciences, 
% Queensland University of Technology, Brisbane, Australia.

% Reference: Milton Abramowitz and Irene A Stegun (1965) Handbook of 
% Mathematical Functions: With Formulas, Graphs, and Mathematical Tables.
% Dover Publications, New York, 1965.

% This function gives the derivative Weierstrass elliptic P function. 
% The calculation is done using the relations between the Weierstrass
% elliptic function P and the the Jacobi's elliptic functions.
% The input z is one complex number or a vector of complex numbers.
% The input g2 and g3 are the invariants.
% The output P is one complex number or a vector of complex numbers 
% depending of the size of z.
function Pprime = WeierstrassPPrime(z,g2,g3)

    % Using homogeneity relation in Equation (18.2.12) when g3 <  0 
    gsign = 1;
    if (g3<0)
        g3=-1*g3;
        gsign = (-1*1i);
        z = 1i*z;
    end
    
    % Finding the roots e_i of the equation 4*e_i^3-g2*e_i-g3=0
    r = roots([4 0 -g2 -g3]);
    e1 = r(1);
    e2 = r(2);
    e3 = r(3);
    
    % Finding the determinant in Equation (18.1.8)
    determinant = g2^3-27*g3^2;
        
    % If the determinant is positive
    if determinant > 0
        % Calculating the parameter m as of the Jacobi's elliptic function
        % in Equation (18.9.9)
        m = (e2-e3)/(e1-e3);        
        % Calculating the input z* of the Jacobi's elliptic function
        % using Equation (18.9.12)
        zstar = sqrt(e1-e3)*z;
        % Calculating the derivative of the Weierstrass elliptic P function
        % using the Jacobi's elliptic functions and Equation (18.9.12)
        Pprime = -2*(e1-e3)^(3/2)*jacobiCN(zstar,m).*jacobiDN(zstar,m)...
                ./(jacobiSN(zstar,m)).^3;
    end
    
    % If the determinant is negative
    if determinant < 0
        % Calculating H2 in Equation (18.3.5)
        H2 = sqrt((e2-e3)*(e2-e1));
        % Calculating the parameter m as of the Jacobi's elliptic function
        % in Equation (18.9.9)
        m = 1/2-3*e2/(4*H2);
        % Calculating the input z* of the Jacobi's elliptic function
        % using Equation (18.9.12)
        zprime = 2*z*sqrt(H2);
        % Calculating the derivative of the Weierstrass elliptic P function
        % using the Jacobi's elliptic functions and Equation (18.9.12)
        Pprime = -4*H2^(3/2)*jacobiSN(zprime,m).*jacobiDN(zprime,m)...
            ./(1-jacobiCN(zprime,m)).^2;
    end
    % Multiplying P by -i if g3<0 as in the homogeneity relation in
    % Equation (18.2.12)
     Pprime = Pprime*gsign;
end