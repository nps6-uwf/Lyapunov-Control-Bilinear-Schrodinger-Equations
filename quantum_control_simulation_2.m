% Nick Sebasco
% May 1st 2023

clc; clear;

% Differential Schrondinger equation:
% p = psi, [x]dot = dx/dt for some x, then,
% p[dot] = -ip(H0 + uH1 + omega)

% n-state quantum system
n = 3;

% states
syms p1 p2 p3;
psi = [p1 p2 p3].'


% Hermitian matrices
H0 = [0 0 0; 0 1 0; 0 0 3/2];
H1 = [0 1 0; 1 0 1; 0 1 0];

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

% choose an eigenvector/ eigenvalue state of interest
idx = 1;
lambda = D(idx,idx)
phi = V(:,idx)

% fictious control
omega = -lambda - b * imag(dot(psi, phi))
u = -a*imag(dot(H1*psi, phi))

pdot = -1i*(H0 + u*H1 + omega * eye(3)) * psi

syms pf1(t) pf2(t) pf3(t);

%pdot = subs(pdot, p1, pf1)
%pdot = subs(pdot, p2, pf3)
%pdot = subs(pdot, p3, pf3)

% time interval of interest
tspan = [0;45];

func = @(x,y,z) subs([pdot(1,:); pdot(2,:); pdot(3,:)],{p1,p2,p3},{x,y,z});

% this is the correct anonymous function that will be supplied to my
% numerical ode solver.  I just copied the rows from pdot and replaced
% {p1,p2,p3} with {x,y,z}, for some reason using the subs function did not
% work so I did this manually.
func_3 = @(x,y,z) [
             - y*((imag(y)*1i)/2 + (imag(z)*1i)/2) - z*((imag(y)*1i)/2 + (imag(z)*1i)/2) - (x*imag(x)*1i)/2;
  - y*((imag(x)*1i)/2 + 1i) - x*((imag(y)*1i)/2 + (imag(z)*1i)/2) - z*((imag(y)*1i)/2 + (imag(z)*1i)/2);
- z*((imag(x)*1i)/2 + 3i/2) - x*((imag(y)*1i)/2 + (imag(z)*1i)/2) - y*((imag(y)*1i)/2 + (imag(z)*1i)/2)
];

func_2 = @(x,y,z) [
                            - (x*imag(x)*1i)/2 - (y*imag(y)*1i)/2;
- y*((imag(x)*1i)/2 + 1i) - (x*imag(y)*1i)/2 - (z*imag(y)*1i)/2;
                   - z*((imag(x)*1i)/2 + 3i/2) - (y*imag(y)*1i)/2    
];

I0 = [0; 1/sqrt(2); 1/sqrt(2)];


I00 = [1; 20; 90];
I0 = 1/norm(I00) * I00;

func(1,1,1);
func_2(1,1,1);

% numerical ode solver ode15s
% https://www.mathworks.com/matlabcentral/answers/1672109-unable-to-find-symbolic-solution-warning
[T,SOL] = ode15s(@(t,sol)func_2(sol(1),sol(2),sol(3)),tspan,I0);

psi_1 = SOL(:,1)

% I need to compute (elemnet-wise) complex magnitude squared, kind of
% tedious, MATLAB should have a built-in for this.
for ii=1: length(psi_1)
psi_1(ii)=abs(psi_1(ii))^2;
end

ustar = subs(u, {p2,p3}, {SOL(:,2), SOL(:,3)});
plot(T,psi_1,'Color',[0.5 0.6 0],'LineWidth',2)
hold on
plot(T,ustar,'--','Color',[0.5 0 0.8],'LineWidth',2)
hold on

% not shown in actual paper (convergence of fictious control)
omegastar = subs(omega, p1, SOL(:,1));
plot(T, omegastar,'--','Color',[0.7 0.2 0],'LineWidth',1)

grid on
legend('$|\Psi_1|^2$','u','$\omega$','Interpreter','latex')
xlabel('Time (s)') 

title(sprintf('a=%.1f, b=%.1f',a,b))




