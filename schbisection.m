clc
clear 
dx=0.000005;
N=10000000;
x=-N*dx/2+dx:dx:N*dx/2;
U=ones(1,N);
S1=N/2-100000; S2=N/2+100000; L=S2-S1; dx=1/L;
gamma=1/3;
g_min=0; g_step=0.01; g_max=0.13;
E=zeros(1,length(g_min:g_step:g_max);
n=0;
for gamma=g_min:g_step:g_mmax
    n=n+1;
    gamma
U(S1:S2)=1;
U(S1:S2)=cantor(U(S1:S2),14,0,gamma);
%U(S1:S2)=0;
%plot(x,U)
%axis([-0.5 0.5 -0.1 1.1]);
%------------------------------------------------------------------------%
it=0;
E1=0.8;
E2=1;
W=zeros(1,N);
K=zeros(1,N);
% updating the wave function---------------------------------------------%
while 1
    W(1)=0; % initial condition
    W(2)=1;
    K(1:N)=(E1-U(1:N));
    for i=2:N-1
        W(i+1)=(2*(1-5/12*dx*dx*K(i))*W(i)-(1+1/12*dx*dx*K(i-1))*W(i-1))/(1+1/12*dx*dx*K(i+1));
    end
    W1=W(N);
    W(1)=0; % initial condition
    W(2)=1;
    K(1:N)=(E2-U(1:N));
    for i=2:N-1
        W(i+1)=(2*(1-5/12*dx*dx*K(i))*W(i)-(1+1/12*dx*dx*K(i-1))*W(i-1))/(1+1/12*dx*dx*K(i+1));
    end
    W2=W(N);
    W(1)=0; % initial condition
    W(2)=1;
    K(1:N)=((E1+E2)/2-U(1:N));
    for i=2:N-1
        W(i+1)=(2*(1-5/12*dx*dx*K(i))*W(i)-(1+1/12*dx*dx*K(i-1))*W(i-1))/(1+1/12*dx*dx*K(i+1));
    end
    W3=W(N);
    if W3*W1<0
        E2=(E1+E2)/2;
    end
    if W3*W2<0
        E1=(E1+E2)/2;
    end
    if abs(E1-E2)<1e-10
        break;
    end
    %{
    plot(x,W,'LineWidth',2);
    hold on;
    plot(x,U.*max(W));
    hold off;
    axis([-N*dx/2 N*dx/2 min(W) max(W)]);
    grid on
    drawnow;
      %} 
    it=it+1;
    %disp((E1+E2)/2)
end
%E=(E1+E2)/2;
E(n)=(E1+E2)/2;
end
%{
K(1:N)=(E(n)-U(1:N));
for i=2:N-1
    W(i+1)=(2*(1-5/12*dx*dx*K(i))*W(i)-(1+1/12*dx*dx*K(i-1))*W(i-1))/(1+1/12*dx*dx*K(i+1));
end
A=1/trapz(W.*W);
Ws=sqrt(A).*W;
format long
disp(it)
% Normalizing probability density funtion of the particle----------------%
plot(x,Ws,'LineWidth',2)
grid on;
%}