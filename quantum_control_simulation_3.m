% Nick Sebasco
% May 1st 2023

clc; clear;

% Differential Schrondinger equation:
% p = psi, [x]dot = dx/dt for some x, then,
% p[dot] = -ip(H0 + uH1 + omega)

% n-state quantum system
n = 3;

% States
psi = sym('p%d',[1,n]).'
syms(['pf(t)'],[n,1])


% Hermitian matrices
H0 = randomComplexHermitian(n)
H1 = randomComplexHermitian(n)

% Verify Hermiticity
% fprintf("H0 Hermitian: %s\n",mat2str(all(all(ctranspose(H0) == H0))))
% fprintf("H1 Hermitian: %s\n",mat2str(all(all(ctranspose(H1) == H1))))
if all(all(ctranspose(H0) == H0))
   fprintf("Hermiticity has been verified for H0 and H1.\n") 
end

% Find the eigenvalues and corresponding eigenvectors of H0.
% D is a diagonal matrix of eigenvalues, V is a matric such
% that the columns are the eigenvectors.
[V, D] = eig(H0);

% Enumerate eigenvalues and corresponding eigenvectors.
% Verify that the eigenvectors are normalized.  If they
% are not normalize them.
for idx=1:n
    D(idx,idx);
    V(:,idx);
    if norm(V(:,idx)) == 1
        continue
    else
        fprintf("eigenvector (%d) undergoing normalization...\n", idx);
        V(:,idx) = 1/norm(V(:,idx)) * V(:,idx);
    end
end

fprintf("All eigenvectors have been normalized.\n")

% a & b control parameters must be greater than 0 to satisfy Lyapunov
% requirement.
a = 1/2;
b = 1/2;

% Choose an eigenvector/ eigenvalue state of interest
idx = 1;
lambda = D(idx,idx)
phi = V(:,idx)

% Fictious control
omega = -lambda - b * imag(dot(psi, phi))
u = -a*imag(dot(H1*psi, phi))

pdot = -1i*(H0 + u*H1 + omega * eye(3)) * psi

% Time interval of interest
tspan = [0;45];

% A more robust way of transforming pdot into a function.
func = matlabFunction(pdot, "Vars", psi);

% Initial condition
I00 = [0; 1/sqrt(2); 1/sqrt(2)];
I0 = 1/norm(I00) * I00;

% Numerical ode solver ode15s
[T,SOL] = ode15s(@(t,sol)func(sol(1),sol(2),sol(3)),tspan,I0);

psi_1 = SOL(:,1)

% I need to compute (elemnet-wise) complex magnitude squared, kind of
% tedious, MATLAB should have a built-in for this.
for ii=1: length(psi_1)
   psi_1(ii)=abs(psi_1(ii))^2;
end

ustar = subs(u, {psi(1), psi(2),psi(3)}, {SOL(:,1), SOL(:,2), SOL(:,3)});
plot(T,psi_1,'Color',[0.5 0.6 0],'LineWidth',2)
hold on

plot(T,ustar,'--','Color',[0.5 0 0.8],'LineWidth',2)
hold on

% not shown in actual paper (convergence of fictious control)
omegastar = subs(omega, {psi(1), psi(2), psi(3)}, {SOL(:,1), SOL(:,2), SOL(:,3)});
plot(T, omegastar,'--','Color',[0.7 0.2 0],'LineWidth',1)

grid on
legend('$|\Psi_1|^2$','u','$\omega$','Interpreter','latex')
xlabel('Time (s)') 

title(sprintf('a=%.1f, b=%.1f',a,b))
