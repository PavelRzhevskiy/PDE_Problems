%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Forward Euler scheme
%%%%%% Ghost-point approach for Neumann BC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all


%Time and space steps:
dt = 1.e-4; dx = 0.02;
mu = dt/dx^2;

%Boundaries and # of pts:
x0 = 0; xJ = 1;
t0 = 0; tN = 5;
N = (tN - t0) / dt + 1; 
J = (xJ - x0) / dx + 3; %+1 ghost point at each boundary

%initial conditions
%Note that x(1) - left ghost point, x(2) - left boundary
%x(J-1) - right boundary, x(J) - right ghost point
x = x0-dx : dx : xJ+dx;
f = (sin(pi*x)).^2;
u = zeros(J, N);
u(:,1) = f; 
u(1,1) = u(3,1); u(J-2,1) = u(J,1); % Introduce ghost points at t = 0

%Loop to find approximate solution:
for n=1:N-1
    %Internal pts (j = 2..N-1)
    for j=2:J-1
        u(j,n+1) = u(j,n) + mu*(u(j+1, n) - 2 * u(j, n) + u(j-1, n)) + ...
            dt * u(j,n) * (1 - u(j,n));
    end 
    %Ghost pts (j=1,J):
    u(1,n+1) = u(3,n+1);
    u(J,n+1) = u(J-2, n+1);
end

%Plot results
tt = 0 : dt : (N-1)*dt;
fig1 = figure(1);

colormap("winter")
surf(tt,x(2:J-1),u(2:J-1,:),'EdgeColor','none');
xlabel('t')
ylabel('x')
zlabel('u')
xlim([0,5])
colorbar
fontsize(fig1, 15, "points")

fig2 = figure(2)
uSlice = -log(1-u(round(J/2),10000:N));
plot(tt(10000:N),uSlice, 'Color','blue', 'LineWidth', 2);
xlim([1,5])
xlabel('t');
ylabel('-log(1-U(0.5, t))');

%Fitting data with pol1
fontsize(fig2, 16, "points")
c=polyfit(tt(10000:N), uSlice, 1)
y_est = polyval(c,tt(10000:N));
hold on
plot(tt(10000:N),y_est,'r--','LineWidth',3)
hold off




