clear
clc
close all

N=41; %number of points
length=3; %length of nozzle
C=0.5;  
dx=length/(N-1); 
cx=0.01;


rho=zeros(1,N); 
V=zeros(1,N);
T=zeros(1,N);

x=linspace(0, length, N); %non-dimensionalized
A=(1+2.2*(x-1.5).^2)/1;

%intial conditions t=0
rho((-1<x)&(x<0.5))=1;
rho((0.4999999999<x)&(x<1.5))=1-0.366*(x((0.4999999999<x)&(x<1.5))-0.5);
rho((1.499999<x)&(x<2.1))=0.634-0.702*(x((1.499999<x)&(x<2.1))-1.5);
rho((2.0999<x)&(x<3.1))=0.5892-0.10228*(x((2.0999<x)&(x<3.1))-2.1);


T((-1<x)&(x<0.5))=1;
T((0.4999999999<x)&(x<1.5))=1-0.167*(x((0.4999999999<x)&(x<1.5))-0.5);
T((1.499999<x)&(x<2.1))=0.833-0.4908*(x((1.499999<x)&(x<2.1))-1.5);
T((2.0999<x)&(x<3.1))=0.93968-0.0622*(x((2.0999<x)&(x<3.1))-2.1);

V=(0.59)./(rho.*A);
total_steps=15000;

p = rho.*T;

g=1.4;
R=287; %J/KGK, exoume aera.

rho_15=zeros(1,total_steps);
V_15=zeros(1,total_steps);
T_15=zeros(1,total_steps);
P_15=zeros(1,total_steps);

drdt_avg_16=zeros(1,total_steps);
dVdt_avg_16=zeros(1,total_steps);
dTdt_avg_16=zeros(1,total_steps);


%V1 changes with time, it floats!!!
%Vend, rho.end, Tend also float!!!
step=1;

    i_shock=round(2.1/dx);

while step < 1000
    step
    
    a = sqrt(T);
    dt = min(C * dx ./ (a + V))

    if imag(dt) ~= 0 | dt<0
        disp('Αρνητικος ή φανταστικος χρονος')
        break
    end
    
    rho_old = rho;
    V_old = V;
    T_old = T;

    % Predictor step for i=1:i_shock-1
    for i = 2:i_shock-1
        drdt_pre(i) = -rho(i) * (V(i+1) - V(i)) / dx - rho(i) * V(i) * (log(A(i+1)) - log(A(i))) / dx ...
                      - V(i) * (rho(i+1) - rho(i)) / dx;
        dVdt_pre(i) = -V(i) * (V(i+1) - V(i)) / dx - 1/g * ((T(i+1) - T(i)) / dx + T(i) / (rho(i)) * (rho(i+1) - rho(i)) / dx);
        dTdt_pre(i) = -V(i) * (T(i+1) - T(i)) / dx - (g - 1) * T(i) * ((V(i+1) - V(i)) / dx + V(i) * (log(A(i+1)) - log(A(i))) / dx);
        S1(i) = cx * abs(p(i+1) - 2 * p(i) + p(i-1)) / (p(i+1) + 2 * p(i) + p(i-1)) * (rho(i+1) - 2 * rho(i) + rho(i-1));
        S2(i) = cx * abs(p(i+1) - 2 * p(i) + p(i-1)) / (p(i+1) + 2 * p(i) + p(i-1)) * (V(i+1) - 2 * V(i) + V(i-1));
        S3(i) = cx * abs(p(i+1) - 2 * p(i) + p(i-1)) / (p(i+1) + 2 * p(i) + p(i-1)) * (T(i+1) - 2 * T(i) + T(i-1));
    end

    for i = 2:i_shock-1
        rho(i) = rho(i) + drdt_pre(i) * dt + S1(i);
        V(i) = V(i) + dVdt_pre(i) * dt + S2(i);
        T(i) = T(i) + dTdt_pre(i) * dt + S3(i);
    end

        %corrector
    for i=2:i_shock-1
        
        drdt_corr(i)=-rho(i)*(V(i)-V(i-1))/dx-rho(i)*V(i)*(log(A(i))-log(A(i-1)))/dx...
                  -V(i)*(rho(i)-rho(i-1))/dx;

        dVdt_corr(i)=-V(i)*(V(i)-V(i-1))/dx-1/g*((T(i)-T(i-1))/dx...
                  +T(i)/(rho(i))*(rho(i)-rho(i-1))/dx);

        dTdt_corr(i)=-V(i)*(T(i)-T(i-1))/dx-(g-1)*T(i)*((V(i)-V(i-1))/dx...
                  +V(i)*(log(A(i))-log(A(i-1)))/dx);


        S1(i)=cx*abs(p(i+1)-2*p(i)+p(i-1))/(p(i+1)+2*p(i)+p(i-1))*(rho(i+1)-2*rho(i)+rho(i-1));
        S2(i)=cx*abs(p(i+1)-2*p(i)+p(i-1))/(p(i+1)+2*p(i)+p(i-1))*(V(i+1)-2*V(i)+V(i-1));
        S3(i)=cx*abs(p(i+1)-2*p(i)+p(i-1))/(p(i+1)+2*p(i)+p(i-1))*(T(i+1)-2*T(i)+T(i-1));
    end

    drdt_avg=0.5*(drdt_pre+drdt_corr);
    dVdt_avg=0.5*(dVdt_pre+dVdt_corr);
    dTdt_avg=0.5*(dTdt_pre+dTdt_corr);

    for i=2:i_shock-1

        rho(i)=rho_old(i)+drdt_avg(i)*dt+S1(i);
        V(i)=V_old(i)+dVdt_avg(i)*dt+S2(i);
        T(i)=T_old(i)+dTdt_avg(i)*dt+S3(i);

    end


    %kroustiko kima sinthikes
    M(i_shock) = 2.07;
    p(i_shock)=(1+(g-1)/2*M(i_shock)^2)^(-g/(g-1));
    rho(i_shock)=(1+(g-1)/2*M(i_shock)^2)^(-1/(g-1));
    T(i_shock)=(1+(g-1)/2*M(i_shock)^2)^(-1);


    M(i_shock+1) = sqrt((1 + ((g - 1) / 2) * M(i_shock)^2) / (g * M(i_shock)^2 - (g - 1) / 2));
    rho(i_shock+1)=((g+1)*M(i_shock)^2)/(2+(g-1)*M(i_shock)^2)*rho(i_shock);
    p(i_shock+1)=(1+(2*g/(g+1))*(M(i_shock)^2-1))*p(i_shock);
    T(i_shock+1)=((1+(2*g/(g+1))*(M(i_shock)^2-1))*((g+1)*M(i_shock)^2)/(2+(g-1)*M(i_shock)^2))*T(i_shock);
    V(i_shock+1)=M(i_shock+1)*T(i_shock)^0.5;


        % Predictor step for i=i_shock:N-1
    for i = i_shock+2:N-1
        drdt_pre(i) = -rho(i) * (V(i+1) - V(i)) / dx - rho(i) * V(i) * (log(A(i+1)) - log(A(i))) / dx ...
                      - V(i) * (rho(i+1) - rho(i)) / dx;
        dVdt_pre(i) = -V(i) * (V(i+1) - V(i)) / dx - 1/g * ((T(i+1) - T(i)) / dx + T(i) / (rho(i)) * (rho(i+1) - rho(i)) / dx);
        dTdt_pre(i) = -V(i) * (T(i+1) - T(i)) / dx - (g - 1) * T(i) * ((V(i+1) - V(i)) / dx + V(i) * (log(A(i+1)) - log(A(i))) / dx);
        S1(i) = cx * abs(p(i+1) - 2 * p(i) + p(i-1)) / (p(i+1) + 2 * p(i) + p(i-1)) * (rho(i+1) - 2 * rho(i) + rho(i-1));
        S2(i) = cx * abs(p(i+1) - 2 * p(i) + p(i-1)) / (p(i+1) + 2 * p(i) + p(i-1)) * (V(i+1) - 2 * V(i) + V(i-1));
        S3(i) = cx * abs(p(i+1) - 2 * p(i) + p(i-1)) / (p(i+1) + 2 * p(i) + p(i-1)) * (T(i+1) - 2 * T(i) + T(i-1));
    end

    %corrector
    for i=i_shock+2:N-1
        
        drdt_corr(i)=-rho(i)*(V(i)-V(i-1))/dx-rho(i)*V(i)*(log(A(i))-log(A(i-1)))/dx...
                  -V(i)*(rho(i)-rho(i-1))/dx;

        dVdt_corr(i)=-V(i)*(V(i)-V(i-1))/dx-1/g*((T(i)-T(i-1))/dx...
                  +T(i)/(rho(i))*(rho(i)-rho(i-1))/dx);

        dTdt_corr(i)=-V(i)*(T(i)-T(i-1))/dx-(g-1)*T(i)*((V(i)-V(i-1))/dx...
                  +V(i)*(log(A(i))-log(A(i-1)))/dx);


        S1(i)=cx*abs(p(i+1)-2*p(i)+p(i-1))/(p(i+1)+2*p(i)+p(i-1))*(rho(i+1)-2*rho(i)+rho(i-1));
        S2(i)=cx*abs(p(i+1)-2*p(i)+p(i-1))/(p(i+1)+2*p(i)+p(i-1))*(V(i+1)-2*V(i)+V(i-1));
        S3(i)=cx*abs(p(i+1)-2*p(i)+p(i-1))/(p(i+1)+2*p(i)+p(i-1))*(T(i+1)-2*T(i)+T(i-1));
    end

    drdt_avg=0.5*(drdt_pre+drdt_corr);
    dVdt_avg=0.5*(dVdt_pre+dVdt_corr);
    dTdt_avg=0.5*(dTdt_pre+dTdt_corr);

    for i=i_shock+2:N-1

        rho(i)=rho_old(i)+drdt_avg(i)*dt+S1(i);
        V(i)=V_old(i)+dVdt_avg(i)*dt+S2(i);
        T(i)=T_old(i)+dTdt_avg(i)*dt+S3(i);

    end

    %BOUNDARY CONDITIONS
    %EISODOS
    V(1)=2*V(2)-V(3); %floats

    %EKSODOS
    V(end)=2*V(end-1)-V(end-2);         %floats
    T(end)=2*T(end-1)-T(end-2);         %floats
    rho(end)=0.6784/T(end);   %floatsfind f



    step = step + 1;
end


figure(1)
plot(x,rho)
title('Density')
grid on
xlabel('X')