%% 1D wave-propagation solver
function [time, u_record] = wavePropagation1d(input)
    L = input.L; Nx = input.Nx; T = input.T; f0 = input.f0;
    rho = input.rho; c0 = input.c0; c_d = input.c_d;
    w = input.damage_width; Dv = input.D_value;
    dc = input.damage_center;
    x = linspace(0,L,Nx)'; dx = x(2)-x(1);
    D = zeros(Nx,1);
    idx = (x >= dc-w/2)&(x<=dc+w/2);
    D(idx)=Dv;
    E = rho*c0^2; E_eff=(1-D)*E;
    c = sqrt(E_eff/rho);
    dt = 0.9*dx/max(c); Nt = floor(T/dt); dt=T/Nt;
    u_old=zeros(Nx,1); u=zeros(Nx,1);
    ricker=@(t)(1-2*(pi*f0*(t-1/f0)).^2).*exp(-(pi*f0*(t-1/f0)).^2);
    E_face=(E_eff(1:end-1)+E_eff(2:end))/2;
    u_record = zeros(Nx,Nt);
    for n=1:Nt
        t=n*dt; u(1)=ricker(t);
        u_new=zeros(Nx,1);
        for i=2:Nx-1
            fr=E_face(i)*(u(i+1)-u(i))/dx;
            fl=E_face(i-1)*(u(i)-u(i-1))/dx;
            Lu=(fr-fl)/dx;
            damp=c_d*(u(i)-u_old(i))/dt;
            u_new(i)=2*u(i)-u_old(i)+(dt^2/rho)*(Lu-damp);
        end
        u_new(end)=u_new(end-1);
        u_old=u; u=u_new;
        u_record(:,n)=u;
    end
    time = (1:Nt)*dt;
end