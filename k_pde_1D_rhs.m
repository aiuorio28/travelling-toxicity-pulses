function f = k_pde_1D_rhs(t,y,par)
% returns RHS for Klausmeier autoxicity system:
% u_t = u_xx+A(1-u)-uv^2
% v_t = eps^2v_xx+uv^2(1-kv)-Bv-Hvs
% s_t = 1/D(-s+Bv+Hvs)

% parameters
N=par.N;
A = par.A;
B = par.B;
D = par.D;
H = par.H;
eps = par.eps;
k = par.k;

% solution
U = y(1:N);
V = y(N+1:2*N);
S = y(2*N+1:3*N);

% differentiation matrices
Dx = par.Dx;
D2x = par.D2x;

%% RHS KA
f1 = D2x*U+A*(1-U)-U.*V.^2;
f2 = eps^2*D2x*V+U.*V.^2.*(1-k*V)-B*V-H*V.*S;
f3 = 1/D*(-S+B*V+H*V.*S);

f = [f1; f2; f3];


 