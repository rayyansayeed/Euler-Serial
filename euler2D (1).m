function euler2D(N,M)
% euler2D(N,M)
% simulate linearized Euler equations
% N = number of grid points in both dimensions
% M = number of time steps

gam = 1.4;
T = 2;
dt = T/M;
x=linspace(-1,1,N+1)';
x=x(1:end-1);
[X,Y] = meshgrid(x,x);
X=X'; Y=Y';
dx = x(2)-x(1);
rho = 2/gam*exp(-100*(X.^2+Y.^2));
u = zeros(size(X));
v = zeros(size(X));
p = 2*exp(-100*(X.^2+Y.^2));
V = [rho;u;v;p];

MLx11 = eye(N,N);
MLx12 = diag(ones(N-1,1),1)*dt/4/dx - diag(ones(N-1,1),-1)*dt/4/dx;
MLx12(1,N) = -dt/4/dx;
MLx12(N,1) = dt/4/dx;
MLx13 = zeros(N,N);
MLx14 = zeros(N,N);
MLx21 = zeros(N,N);
MLx22 = eye(N,N);
MLx23 = zeros(N,N);
MLx24 = MLx12;
MLx31 = zeros(N,N);
MLx32 = zeros(N,N);
MLx33 = eye(N,N);
MLx34 = zeros(N,N);
MLx41 = zeros(N,N);
MLx42 = gam*MLx12;
MLx43 = zeros(N,N);
MLx44 = eye(N,N);
MLx = [MLx11,MLx12,MLx13,MLx14;MLx21,MLx22,MLx23,MLx24; ...
    MLx31,MLx32,MLx33,MLx34;MLx41,MLx42,MLx43,MLx44];

MRx11 = eye(N,N);
MRx12 = -diag(ones(N-1,1),1)*dt/4/dx + diag(ones(N-1,1),-1)*dt/4/dx;
MRx12(1,N) = dt/4/dx;
MRx12(N,1) = -dt/4/dx;
MRx13 = zeros(N,N);
MRx14 = zeros(N,N);
MRx21 = zeros(N,N);
MRx22 = eye(N,N);
MRx23 = zeros(N,N);
MRx24 = MRx12;
MRx31 = zeros(N,N);
MRx32 = zeros(N,N);
MRx33 = eye(N,N);
MRx34 = zeros(N,N);
MRx41 = zeros(N,N);
MRx42 = gam*MRx12;
MRx43 = zeros(N,N);
MRx44 = eye(N,N);
MRx = [MRx11,MRx12,MRx13,MRx14;MRx21,MRx22,MRx23,MRx24; ...
    MRx31,MRx32,MRx33,MRx34;MRx41,MRx42,MRx43,MRx44];

MLy11 = eye(N,N);
MLy12 = zeros(N,N);
MLy13 = diag(ones(N-1,1),1)*dt/4/dx - diag(ones(N-1,1),-1)*dt/4/dx;
MLy13(1,N) = -dt/4/dx;
MLy13(N,1) = dt/4/dx;
MLy14 = zeros(N,N);
MLy21 = zeros(N,N);
MLy22 = eye(N,N);
MLy23 = zeros(N,N);
MLy24 = zeros(N,N);
MLy31 = zeros(N,N);
MLy32 = zeros(N,N);
MLy33 = eye(N,N);
MLy34 = MLy13;
MLy41 = zeros(N,N);
MLy42 = zeros(N,N);
MLy43 = gam*MLy13;
MLy44 = eye(N,N);
MLy = [MLy11,MLy12,MLy13,MLy14;MLy21,MLy22,MLy23,MLy24; ...
    MLy31,MLy32,MLy33,MLy34;MLy41,MLy42,MLy43,MLy44];

MRy11 = eye(N,N);
MRy12 = zeros(N,N);
MRy13 = -diag(ones(N-1,1),1)*dt/4/dx + diag(ones(N-1,1),-1)*dt/4/dx;
MRy13(1,N) = dt/4/dx;
MRy13(N,1) = -dt/4/dx;
MRy14 = zeros(N,N);
MRy21 = zeros(N,N);
MRy22 = eye(N,N);
MRy23 = zeros(N,N);
MRy24 = zeros(N,N);
MRy31 = zeros(N,N);
MRy32 = zeros(N,N);
MRy33 = eye(N,N);
MRy34 = MRy13;
MRy41 = zeros(N,N);
MRy42 = zeros(N,N);
MRy43 = gam*MRy13;
MRy44 = eye(N,N);
MRy = [MRy11,MRy12,MRy13,MRy14;MRy21,MRy22,MRy23,MRy24; ...
    MRy31,MRy32,MRy33,MRy34;MRy41,MRy42,MRy43,MRy44];

figure(1)
subplot(4,1,1);
surf(X,Y,V(1:N,1:N))
axis([-1,1,-1,1,-3,3]);
subplot(4,1,2);
surf(X,Y,V(N+1:2*N,1:N))
axis([-1,1,-1,1,-2,2]);
subplot(4,1,3);
surf(X,Y,V(2*N+1:3*N,1:N));
axis([-1,1,-1,1,-2,2]);
subplot(4,1,4);
surf(X,Y,V(3*N+1:4*N,1:N))
axis([-1,1,-1,1,-2,2]);
drawnow

for t=dt:dt:T

    % forward in x
    V = MRx*V;
    
    % backward in y
    W = [V(1:N,:)';V(N+1:2*N,:)';V(2*N+1:3*N,:)';V(3*N+1:4*N,:)'];
    W = MLy\W;
    
    % forward in y
    W = MRy*W;
    
    % backward in x
    V = [W(1:N,:)';W(N+1:2*N,:)';W(2*N+1:3*N,:)';W(3*N+1:4*N,:)'];
    V = MLx\V;
    
    subplot(2,2,1);
    surf(X,Y,V(1:N,1:N))
    title(sprintf('Time = %f',t));
    axis([-1,1,-1,1,-1,1]);
    axis equal
    title('\rho')
    subplot(2,2,2);
    surf(X,Y,V(N+1:2*N,1:N))
    axis([-1,1,-1,1,-2,2]);
    axis equal
    title('u')
    subplot(2,2,3);
    surf(X,Y,V(2*N+1:3*N,1:N));
    axis([-1,1,-1,1,-2,2]);
    axis equal
    title('v')
    subplot(2,2,4);
    surf(X,Y,V(3*N+1:4*N,1:N))
    axis([-1,1,-1,1,-2,2]);
    axis equal
    title('p')
    drawnow
    
end

end
    