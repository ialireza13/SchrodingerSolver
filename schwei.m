clc
N=1000000;
U=zeros(1,N);
S1=N/2-10000; S2=N/2+10000; L=S2-S1; dx=1/L;
x=-N*dx/2+dx:dx:N*dx/2;
U(1:S1)=1; U(S2:N)=1;
W=zeros(1,N);
%b=7;
%Energy = zeros(length((1+1.5*pi)/b:(1-(1+1.5*pi)/b)/(b*8):1),2);
jj=0;
prev=0;
hold on
for b=7:2:23
    disp('b=')
    disp(b);
    %for a=0.55:0.55
    for a=(1+1.5*pi)/b:(1-(1+1.5*pi)/b)/(b^1.7):1
        disp('a=')
        disp(a);
        jj=jj+1;
        ii=S1-1;
        for i=-L/2+1:L/2
            ii=ii+1;
            Wei=0;
            for n=0:100
                Wei=Wei+ (a.^n * cos(2*b.^n * pi * i*dx));
            end
            U(ii)=Wei;
        end
        U(S1:S2) = (U(S1:S2)-min(U(S1:S2)))/(max(U(S1:S2))-min(U(S1:S2)));
        %U(S1:S2) = U(S1:S2)./U(N/2);
        %plot(x,U)
        %axis([-0.5 0.5 0 1])
        %plot(S1*dx:dx:S2*dx,U(S1:S2))
        %plot(U(S1:S2))
        %it=0;
        E1=0.85;
        E2=0.95;
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
            if abs(E1-E2)<1e-9
                break;
            end
            %it=it+1;
        end
        E=(E1+E2)/2;
        K(1:N)=(E-U(1:N));
            for i=2:N-1
                W(i+1)=(2*(1-5/12*dx*dx*K(i))*W(i)-(1+1/12*dx*dx*K(i-1))*W(i-1))/(1+1/12*dx*dx*K(i+1));
            end
        format long
        %A=1/trapz(W.*W);
        %Ws=sqrt(A).*W;
        %subplot(1,2,1);
        %plot(Ws)
        %drawnow;
        [M,I] = max(W(S1:S2));
        if I<9000 || I>11000
            disp('Error!');
            break;
        end
        %title(strcat(strcat('b=',b),strcat('  a=',a)));
        Energy(jj,1)=2+log(a)/log(b); Energy(jj,2)=E;
    end
    %subplot(1,2,2);
    plot(Energy(prev+1:jj,1),Energy(prev+1:jj,2));
    drawnow;
    hold on
    prev=jj;
end