%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Problem 2
%%%%%% CN Scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all

% Declare some parameters
maxIterNum = 10000;

% Time and space steps:
dt = 1.e-4; dx = 1.e-3;
tol = min([dt^2, dx^2]);
mu = dt/dx^2;

% Boundaries and # of pts:
x0 = -1; xJ = 1;
t0 = 0; tN = 0.5;
N = (tN - t0) / dt + 1; 
J = (xJ - x0) / dx + 1;

% initial conditions u(:,1):
x = x0 : dx : xJ;
tt = t0 : dt : tN;
u = zeros(J, N);
for i = 1 : J
    u(i,1) = initVal(x(i));
end

%% CN scheme.
%% We solve tridiagonal system of eq-s
%% using Thomas Algorithm - O(J) operations
%% https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
%% implemented by Tamas Kis https://tamaskis.github.io

% Coefficients of B- and C- matrices:
% Note that number of unknowns - J-2
AVec = -mu * ones(J-3,1);
BVec = (2 + 2*mu) * ones(J-2,1);
CVec = -mu * ones(J-3,1);

CMatrCoeff = [mu, 2 - 2 * mu, mu];
DVec = zeros(J-2,1);

%BC:
u(1,:) = exactSol(x0, tt, tol, maxIterNum);
u(J,:) = exactSol(xJ, tt, tol, maxIterNum);
for n = 2 : N
    uOld = u(2:J-1, n-1);
    t = (n-1)*dt;
    bc_l = mu*(u(1,n-1) + u(1,n));
    bc_r = mu*(u(J,n-1) + u(J,n));
    
    DVec(1) = bc_l + CMatrCoeff(2) * uOld(1) + CMatrCoeff(3) * uOld(2);
    for j = 2 : J-3
        DVec(j) = CMatrCoeff(1) * uOld(j-1) + CMatrCoeff(2) * uOld(j) + ...
            CMatrCoeff(3) * uOld(j+1);
    end
    DVec(J-2) =  CMatrCoeff(1) * uOld(J-3) + CMatrCoeff(2) * uOld(J-2) + bc_r;

    u(2:J-1, n) = tridiagonal_vector(AVec, BVec, CVec, DVec);

end


fig1 = figure(1);
colormap("winter")
surf(tt,x,u,'EdgeColor','none');
xlabel('t')
ylabel('x')
zlabel('uAppr')
colorbar
fontsize(fig1, 15, "points")

uExact = zeros(J,N);
for n = 1 : N
    t = (n-1)*dt;
    uExact(:,n) = exactSol(x, t, tol, maxIterNum);
end

fig2 = figure(2);
colormap("spring")
surf(tt,x,uExact,'EdgeColor','none');
xlabel('t')
ylabel('x')
zlabel('uExact')
colorbar
fontsize(fig2, 15, "points")


fig3 = figure(3);
colormap("spring")
surf(tt,x,uExact - u,'EdgeColor','none');
xlabel('t');
ylabel('x');
zlabel('u');
zlim([0,0.01]);
err_L2 = sqrt(sum( (uExact - u).^2, "all")/(N*J))
err_LInf = max(abs(uExact - u),[],"all")
colorbar
fontsize(fig2, 15, "points")


function u = exactSol(x, t, tol, maxIterNum)
    %Loop to find approximate solution:
    uval = 1/2;
    temp = 10;
    n = 0;
    while (any(abs(temp) > tol) && n < maxIterNum)
        temp = 2*(-1)^n * cos(pi*(2*n+1)*x)/(pi*(2*n+1))*exp(-pi^2*(2*n+1)^2*t);
        uval = uval + temp;
        n = n + 1;
    end
    u=uval;
end

function u = initVal(x)
    eps = 1.e-15;
    if (abs(x) < 1/2)
        u = 1;
    elseif (abs(x-1/2)<eps)
        u = 1/2;
    else 
        u = 0;
    end
end