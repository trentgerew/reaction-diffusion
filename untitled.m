clear
clc

N = 100; % number of grid points
L = 1; % repeat length

x = ([1:N] - ones(1,N)/2) * L / N;
dx = L / N;

k12 = 2;
k21 = 1;

K12 = Kmat(k12,N);
K21 = Kmat(k21,N);

dphi1 = phi1(x(2:N),L) - phi1(x(1:N-1),L);
dphi2 = phi2(x(2:N),L) - phi2(x(1:N-1),L);

F1 = Fmat(dphi1,dx);
F2 = Fmat(dphi2,dx);

B1 = Bmat(dphi1,dx);
B2 = Bmat(dphi2,dx);

L1 = Lmat(F1,B1);
L2 = Lmat(F2,B2);

L1(N,1) = B1(1);
L1(1,N) = F1(end);

L2(N,1) = B2(1);
L2(1,N) = F2(end);

M = [L1 - K12 , K21 ; K12 , L2 - K21];

% Initial Condition
p0 = zeros(2*N,1);
p0(round(N/4)) = 1;

% Solve
tstop = 1;
dt = 1e-4;
t = 1:dt:tstop;

p = zeros(2*N,length(t));
p(:,1) = p0;
for i = 2:length(t)
    p(:,i) = (dt * M + ones(2*N, 2*N)) * p(:,i-1);
end

plot(p(1:N,length(t) ))

%% Functions

% Potentials
function phi = phi1(x,L)
    phi = sin(2 * pi * x / L);
end

function phi = phi2(x,L)
    phi = zeros(1,length(x));
end

% Jump rates
function fmat = Fmat(dphi,dx)
    if dphi == zeros(1,length(dphi))
        fmat = zeros(1,length(dphi));
    else
        fmat = dphi ./ (dx^2 * (exp(dphi) - ones(1,length(dphi))));
    end
end

function bmat = Bmat(dphi,dx)
    if dphi == zeros(1,length(dphi))
        bmat = zeros(1,length(dphi));
    else
        bmat = -dphi ./ (dx^2 * (exp(-dphi) - ones(1,length(dphi))));
    end
end

% Matrix formation
function lmat = Lmat(f,b)
    lmat = diag([-f - b , 0]) + diag(b,-1) + diag(f,1);
end

function kmat = Kmat(k,N)
    kmat = k * eye(N,N);
end