clear
clc
close all

%% Simulation Parameters
N = 100; % number of grid points
tstop = .001; % stopping time

%% Problem Parameters
L = .25; % repeat length

% Rate constants
k12 = 2000;
k21 = 1000;

% Initial Condition
p0 = ones(2*N,1) / N;
%p0(round(N/2)) = 1; % 

% See Functions section below to change potentials

%% Simulation
% Grid
x = ([1:N] - ones(1,N)/2) * L / N; % create grid points
dx = L / N; % grid spacing

% K matrices (diagonal N x N matrices with the values of the rate
% constants along the diagonal) (
K12 = Kmat(k12,N);
K21 = Kmat(k21,N);

% Potentials evaluated along grid (equation 9)
dphi1 = phi1(x(2:N),L) - phi1(x(1:N-1),L);
dphi2 = phi2(x(2:N),L) - phi2(x(1:N-1),L);

% Forward jump rates evaluated at half-grid steps (equation 10)
F1 = Fmat(dphi1,dx);
F2 = Fmat(dphi2,dx);

% Backward jump rates evaluated at half-grid steps (equation 11)
B1 = Bmat(dphi1,dx);
B2 = Bmat(dphi2,dx);

% Form L matrices (equations 13-15)
L1 = Lmat(F1,B1);
L2 = Lmat(F2,B2);

% Specify periodic boundary conditions
L1(N,1) = B1(2); % (equation 24)
L1(1,N) = F1(end-1); % (equation 25)

L2(N,1) = B2(2); % (equation 24)
L2(1,N) = F2(end-1); % (equation 25)

% Form M matrix (equation 23)
M = [L1 - K12 , K21 ; K12 , L2 - K21];

% Solve
tspan = [0 tstop]; % time domain
[t,p] = ode45(@(t,p) odefun(t,p,M), tspan, p0); % solve linear system
p = p'; % tranpose solution so rows are n and cols are t

% Separate solutions
p1 = p(1:N,:);
p2 = p(N+1:end,:);

% Plot results
plt(p1,x,t,1)
plt(p2,x,t,2)

pl(x,p1,1)
pl(x,p2,2)

% Calculate probability over time
rho1 = inte(p1,dx);
rho2 = inte(p2,dx);

% Plot probility over time
figure;
plot(t,rho1,'r-',t,rho2,'b-')
xlabel('t')
legend('\rho_1(t)','\rho_2(t)')

% Show end probability
formatSpec = 'rho_%d (t_f) = %1.8f\n';
fprintf(formatSpec,1,rho1(end))
fprintf(formatSpec,2,rho2(end))

%% Parameter Exports

init = [rho1(1) rho2(1)]; % Initial Conditions

% Export
save('params.mat','init','k12','k21','tstop')

%% Functions

% Potentials
function phi = phi1(x,L) % potential 1
    %phi = sin(2 * pi * x / L) - sin(4 * pi * x / L) / 2;
    phi = zeros(1,length(x));
end

function phi = phi2(x,L) % potential 2
    phi = zeros(1,length(x));
end

% Jump rates
function fmat = Fmat(dphi,dx) % foward rate
    if dphi == zeros(1,length(dphi)) % if the potential is identically 0
        fmat = zeros(1,length(dphi) + 1); % the rate is identically 0
    else
        fmat = [dphi ./ (dx^2 * (exp(dphi) - ones(1,length(dphi)))) , 0]; % (equation 10)
    end
end

function bmat = Bmat(dphi,dx) % backward rate
    if dphi == zeros(1,length(dphi)) % if the potential is identically 0
        bmat = zeros(1,length(dphi) + 1); % the rate is identically 0
    else
        bmat = [0 , -dphi ./ (dx^2 * (exp(-dphi) - ones(1,length(dphi))))]; % (equation 11)
    end
end

% Matrix formation
function lmat = Lmat(f,b) % L matrix
    lmat = diag(-f - b) + diag(b(2:end),1) + diag(f(1:end-1),-1);
end

function kmat = Kmat(k,N) % K matrix
    kmat = k * eye(N,N);
end

% Differential equation
function dpdt = odefun(~,p,M)
    dpdt = M * p;
end

% Result plot
function [] = plt(p,x,t,pltnum) % 3D surface plot p(x,t)
    figure;
    surf(t,x,p)
    title(['Density in space and time, n=' num2str(pltnum)])
    zlabel('p_n(x,t)')
    xlabel('t')
    ylabel('x')
end

function [] = pl(x,p,pltnum) % 2D plot p(x,t_f)
    figure;
    plot(x,p(:,end))
    title(['Density in space at t_f, n=' num2str(pltnum)])
    xlabel('x')
    ylabel('p_n(x,t_f)')
end

% Area
function area = inte(p,dx)
    for i=1:size(p,2)
        area(i) = trapz(p(:,i)) * dx;
    end
end