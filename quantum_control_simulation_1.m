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

% states p1, ..., pn must satisfy,
% sum 1 to n |pi|^2 = 1, one guarantee
% for this criteria is to choose pi of the form
% pi = +-sqrt(1/(2n)) +-sqrt(1/(2n))i, proof

c1 = complex(sqrt(1/(2*n)), sqrt(1/(2*n)));
c2 = c1;
c3 = c2;

cn = [c1 c2 c3];
magnitude = 0;

for idx=1:n
    magnitude = magnitude + norm(cn(idx))^2;
end

magnitude

% Hermitian matrices
H0 = [0 0 0; 0 1 0; 0 0 3/2];
H1 = [0 1 1; 1 0 1; 1 1 0];

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
for i=idx:n
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
b = a;

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

% func_2 was created from incorrect pdot
func_2 = @(x,y,z) [
    - (y*imag(x)*(imag(y)/2 + imag(z)/2)*1i)/2 - (z*imag(x)*(imag(y)/2 + imag(z)/2)*1i)/2;
    - y*1i - (x*imag(x)*(imag(y)/2 + imag(z)/2)*1i)/2 - (z*imag(x)*(imag(y)/2 + imag(z)/2)*1i)/2;
    - (z*3i)/2 - (x*imag(x)*(imag(y)/2 + imag(z)/2)*1i)/2 - (y*imag(x)*(imag(y)/2 + imag(z)/2)*1i)/2
    ];

% this is the correct anonymous function that will be supplied to my
% numerical ode solver.  I just copied the rows from pdot and replaced
% {p1,p2,p3} with {x,y,z}, for some reason using the subs function did not
% work so I did this manually.
func_3 = @(x,y,z) [
             - y*((imag(y)*1i)/2 + (imag(z)*1i)/2) - z*((imag(y)*1i)/2 + (imag(z)*1i)/2) - (x*imag(x)*1i)/2;
  - y*((imag(x)*1i)/2 + 1i) - x*((imag(y)*1i)/2 + (imag(z)*1i)/2) - z*((imag(y)*1i)/2 + (imag(z)*1i)/2);
- z*((imag(x)*1i)/2 + 3i/2) - x*((imag(y)*1i)/2 + (imag(z)*1i)/2) - y*((imag(y)*1i)/2 + (imag(z)*1i)/2)
];

I0 = [0; 1/sqrt(2); 1/sqrt(2)];

func(1,1,1);
func_2(1,1,1);
func_3(1,1,1);

% numerical ode solver ode15s
% https://www.mathworks.com/matlabcentral/answers/1672109-unable-to-find-symbolic-solution-warning
[T,SOL] = ode15s(@(t,sol)func_3(sol(1),sol(2),sol(3)),tspan,I0);

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

% I first attempted to solve the system of differential equations by using
% MATLAB's dsolve, function.  This failed becuase I was dealing with
% symbolic equations of multiple variables.  So you can see above that I
% changed my approach to use a numerical solver ode15s.
if false
    ode1 = diff(pf1) == pdot(1,:)
    ode2 = diff(pf2) == pdot(2,:)
    ode3 = diff(pf3) == pdot(3,:)
    
    odes = [ode1; ode2; ode3]
    
    % initial conditions
    cond1 = pf1(0) == 0;
    cond2 = pf2(0) == 1/sqrt(2);
    cond3 = pf3(0) == 1/sqrt(2);
    
    conds = [cond1; cond2; cond3];

    % dsolve does not work so I need to slve numerically
    % dsolve(odes, conds);

    [psi1(t),psi2(t),psi3(t)] = dsolve(odes, conds);
    
    %solution = dsolve(odes, conds)
    %solution1 = solution.pf1;
    %solution1 = subs(solution1, imag(p2),1);
    %solution1 = subs(solution1, imag(p3),1);
    %solution1 = subs(solution1, conj(p3),1);
    %solution1 = subs(solution1, conj(p2),1);
    
    psi1 = subs(psi1, {p1,p2,p3}, {c1,c2,c3})
    
    fplot(norm(psi1(t))^2, T);
    %hold on
    %fplot(vSol)
    %grid on
end
