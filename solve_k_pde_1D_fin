function solution = solve_k_pde_1D_fin(tend,K)

% solves Klausmeier PDE with RHS k_pde_1D_rhs.m
% and Jacobian Dk_pde_1D_rhs.m
% plots solution after each time tend for K iterations

par.N = 29970; % (for case i)
% par.N = 10000; % (for case ii)
N = par.N;
dt = tend;

%% setup initial condition

% Case (i) 1: NO plateau 
  % (tend,K)=(100,500)
par.A = 1.5;
par.B = 0.2;
par.D = 3160;
par.H = 0.1;
par.eps = 0.001;
par.k = 1.059;
par.Lx = 10;

% % Case (i) 2: superslow plateau
  %  % (tend,K)=(100,500)
% par.A = 1.5;
% par.B = 0.2;
% par.D = 2277;
% par.H = 0.1;
% par.eps = 0.001;
% par.k = 0.955;
% par.Lx = 20;

% % Case (ii)
  %  % (tend,K)=(100,500)
% par.A = 1.5;
% par.B = 0.2;
% par.D = 37492;
% par.H = 0.1;
% par.eps = 0.01;
% par.k = 0.955;
% par.Lx = 60;

Lx = par.Lx;
par.hx = Lx/(N-1); hx = par.hx;
x = (1:N)'*hx;
par.x = x;

% load('data29970_i_ss_fin') % (initial data for case i with superslow plateau)
load('data29970_i_noss_fin2') % (initial data for case i without superslow plateau)
% load('dataiifin10000') % (initial data for case ii)
% sol = endstate;

sol = Expression1; % (only for case (ii))

%% differentiation matrices

e = ones(N,1);

% d_x
Dx = sparse(1:N-1,[2:N-1 N],ones(N-1,1)/2,N,N);
Dx(1,N) = -1/2; % Periodic boundary conditions
Dx = (Dx - Dx')/hx;
%Dx(1,2) = 0; Dx(N,N-1) = 0; % Neumann BCs

% d_xx
D2x = sparse(1:N-1,[2:N-1 N],ones(N-1,1),N,N) - sparse(1:N,[1:N],e,N,N);
D2x = (D2x + D2x');
%D2x(1,2)=2; D2x(N,N-1)=2; % Neumann Bcs
D2x(1,N) = 1; D2x(N,1) = 1; % Periodic boundary conditions
D2x = D2x/hx^2;

par.Dx = Dx;
par.D2x = D2x;

%% solve PDE

solution= [];
times = [];
Dk_pde_1D_rhs_s = @(t,y)Dk_pde_1D_rhs(t,y,par);

options=odeset('RelTol',1e-8,'AbsTol',1e-8,'Jacobian',Dk_pde_1D_rhs_s);

for j=0:K-1
    j
    time = [0:dt:tend];
    sol = ode15s(@(t,y)k_pde_1D_rhs(t,y,par), time, sol,options);
    times = [times j*tend];
    sol = sol.y(:,end);
    solution = [solution sol];
end

%% save solution

save('ka_pulse_end_casei_noss_fin2_29970','sol'); (final solution for case i without superslow plateau)
% save('ka_pulse_end_casei_ss_fin_29970','sol'); (final solution for case i with superslow plateau)
% save('ka_pulse_end_caseii_fin_10000','sol'); (final solution for case ii)


%% plot solution

% spacetime plot
figure(2)
s = surf(x, times, solution(N+1:2*N,:)');
colormap(flipud(summer));
s.EdgeColor = 'none';
xlim([0 par.Lx])
% zlim([0 3])
%xlabel('x'), ylabel('t'), zlabel('V')

% solution profiles
for i=1:length(times)
    figure(1);
    subplot(1,3,1);
    plot(x,solution(1:N,i));
    title(['U-profile at time t=' num2str(times(i))]);
    subplot(1,3,2);
    plot(x,solution(N+1:2*N,i));
    title(['V-profile at time t=' num2str(times(i))]);
    subplot(1,3,3);
    plot(x,solution(2*N+1:3*N,i));
    title(['S-profile at time t=' num2str(times(i))]);
    pause;
end

% 3d plot in phase space
U = solution(1:N,end);
V = solution(N+1:2*N,end);
S = solution(2*N+1:3*N,end);

figure(3)
plot3(V,Dx*V,S)
xlabel('V')
ylabel('V_x')
zlabel('S')

figure(4)
plot3(V,U,S)
xlabel('V')
ylabel('U')
zlabel('S')
