function Df = Dk_pde_1D_rhs(t,y,par)
% returns Jacobian of RHS for Klausmeier+autotoxicity system

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

Dx = par.Dx;
D2x = par.D2x;

D11 = D2x-A*speye(N,N)-spdiags(V.^2,0,N,N);
D12 = -spdiags(2*U.*V,0,N,N);
D13 = sparse(N,N);
D21 = spdiags(V.^2.*(1-k*V),0,N,N);
D22 = eps^2*D2x+spdiags(U.*(2*V-3*k*V.^2),0,N,N)-B*speye(N,N)-spdiags(H*S,0,N,N);
D23 = -spdiags(H*V,0,N,N);
D31 = sparse(N,N);
D32 = spdiags(B/D*V,0,N,N)+spdiags(H/D*S,0,N,N);
D33 = -1/D*speye(N,N)+spdiags(H/D*V,0,N,N);

Df=[D11,D12, D13;D21,D22, D23;D31,D32, D33];
 end