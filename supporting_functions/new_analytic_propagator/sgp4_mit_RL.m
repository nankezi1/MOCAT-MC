function [sat,r,v,oe_mean_out]=sgp4_mit(sat,time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code implements MIT's version of sgp4 
% Author: Richard Linares, MIT 09/03/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%inputs
no=sat.no;
Omegao=sat.nodeo;
eo=sat.ecco;
io=sat.inclo;
Mo=sat.mo;
omegao=sat.argpo;
Bstar=sat.bstar;

to=0;
t=time;

%constants
radiusearthkm = 6378.135; 
aE=1;%units of 1 Earth Radius 
J3=-0.253881e-5;
mu=398600.79964; 
XKMPER=radiusearthkm;
xke    = 60.0* sqrt(mu/(radiusearthkm*radiusearthkm*radiusearthkm));
vkmpersec     = radiusearthkm * xke/60.0;

% k2=1/2*J2*aE^2;
% k4=-(3/8)*J4*aE^4;
k2=5.413080e-4;
k4=0.62098875e-6;
A30=-J3*aE^3;
qo=1;
ke=xke;%sqrt(mu*(1/XKMPER^2));
s=1.01222928;
qo=(1.88027916e-9)^(1/4)+s;

a1=(ke/no)^(2/3);

delta1=(3/2)*(k2/a1^2)*((3*cos(io)^2-1)/((1-eo^2)^(3/2)));

ao=a1*(1-(1/3)*delta1-delta1^2-(134/81)*delta1^3);

deltao=(3/2)*(k2/ao^2)*((3*cos(io)^2-1)/(1-eo^2)^(3/2));

nopp=no/(1+deltao);

aopp=(ao)/(1-deltao);

if 98<ao*(1-eo)<156
    ss=aopp*(1-eo)-s-aE;
    term=((qo-ss)-(s-ss))^4;
end

if aopp*(1-eo)<98
    ss=20/XKMPER+aE;
    term=((qo-ss)-(s-ss))^4;
end

%This is a term that is needed for power law density model used in SGP$
% term=(qo-s)^4
term=1.88027916e-9;

theta=cos(io);

xi=1/(aopp-s);
betao=sqrt(1-eo^2);
eta=aopp*eo*xi;

C2=term*xi^4*nopp*(1-eta^2)^(-7/2)*(aopp*(1+(3/2)*eta^2+4*eo*eta+eo*eta^3)+...
    +(3/2)*((k2*xi)/(1-eta^2))*(-1/2+3/2*theta^2)*(8+24*eta^2+3*eta^4));

C1=Bstar*C2;

C3=(term*eta^5*A30*nopp*aE*sin(io))/(k2*eo);

C4=2*nopp*term*xi^4*aopp*betao^2*(1-eta^2)^(-7/2)*((2*eta*(1+eo*eta)+1/2*eo+1/2*eta^3)-((2*k2*xi)/(aopp*(1-eta^2)))*...
    (3*(1-3*theta^2)*(1+3/2*eta^2-2*eo*eta-1/2*eo*eta^3)+ 3/4*(1-theta^2)*(2*eta^2-eo*eta-eo*eta^3)*cos(2*omegao)) );

C5=2*term*xi^4*aopp*betao^2*(1-eta^2)^(-7/2)*(1+(11/4)*eta*(eta+eo)+eo*eta^3);

D2=4*aopp*xi*C1^2;

D3=(4/3)*aopp*xi^2*(17*aopp+s)*C1^3;

D4=(2/3)*aopp*xi^3*(221*aopp+31*s)*C1^4;

M_DF=Mo+(1+(3*k2*(-1+3*theta^2))/(2*aopp^2*betao^3)+(3*k2^2*(13-78*theta^2+137*theta^4))/(16*aopp^4*betao^7))*nopp*(t-to);

omega_DF=omegao+(-((3*k2*(1-5*theta^2))/(2*aopp^2*betao^4))+((3*k2^2*(7-114*theta^2+395*theta^4))/(16*aopp^4*betao^8))...
    +((5*k4*(3-36*theta^2+49*theta^4))/(4*aopp^4*betao^8)))*nopp*(t-to);

Omega_DF=Omegao+(-((3*k2*theta)/(aopp^2*betao^4))+...
    ((3*k2^2*(4*theta-19*theta^3))/(2*aopp^4*betao^8))+...
    ((5*k4*theta*(3-7*theta^2))/(2*aopp^4*betao^8)))*nopp*(t-to);

delta_omega=Bstar*C3*(cos(omegao))*(t-to);
delta_M=-(2/3)*term*Bstar*xi^4*((aE)/(eo*eta))*((1+eta*cos(M_DF))^3-(1+eta*cos(Mo))^3);

Mp=M_DF+delta_omega+delta_M;
omega=omega_DF-delta_omega-delta_M;
Omega=Omega_DF-(21/2)*((nopp*k2*theta)/(aopp^2*betao^2))*C1*(t-to)^2;

e=eo-Bstar*C4*(t-to)-Bstar*C5*(sin(Mp)-sin(Mo));

a=aopp*(1-C1*(t-to)-D2*(t-to)^2-D3*(t-to)^3-D4*(t-to)^4)^2;

IL=Mp+omega+Omega+nopp*((3/2)*C1*(t-to)^2+(D2+2*C1^2)*(t-to)^3+...
    1/4*(3*D3+12*C1*D2+10*C1^3)*(t-to)^4+...
    +1/5*(3*D4+12*C1*D3+6*D2^2+30*C1^2*D2+15*C1^4)*(t-to)^5);

beta=sqrt(1-e^2);

n=ke/(a^(3/2));

a_xN=e*cos(omega);

IL_L=(A30*sin(io))/(8*k2*a*beta^2)*(e*cos(omega))*((3+5*theta)/(1+theta));

a_yN_L=(A30*sin(io))/(4*k2*a*beta^2);

IL_T=IL+IL_L;

a_yN=e*sin(omega)+a_yN_L;

U=IL_T-Omega;

termEW=U;
for k=1:15
    DtermEW=(U-a_yN*cos(termEW)+a_xN*sin(termEW)-termEW)/(-a_yN*sin(termEW)-a_xN*cos(termEW)+1);
    termEW=termEW+DtermEW;
end

ecosE=a_xN*cos(termEW)+a_yN*sin(termEW);
esinE=a_xN*sin(termEW)-a_yN*cos(termEW);

e_L=sqrt(a_xN^2+a_yN^2);
p_L=a*(1-e_L^2);

r=a*(1-ecosE);

r_dot=ke*(sqrt(a)/r)*esinE;
rf_dot=ke*(sqrt(p_L)/r);

cosu=(a/r)*(cos(termEW)-a_xN+(a_yN*(esinE))/(1+sqrt(1-e_L^2)));

sinu=(a/r)*(sin(termEW)-a_yN-(a_xN*(esinE))/(1+sqrt(1-e_L^2)));

u=atan2(sinu,cosu);

Dr=(k2/(2*p_L))*(1-theta^2)*cos(2*u);

Du=-(k2/(2*p_L^2))*(7*theta^2-1)*sin(2*u);

DOmega=((3*k2*theta)/(2*p_L^2))*sin(2*u);

Di=((3*k2*theta)/(2*p_L^2))*cos(2*u)*sin(io);

Dr_dot=-((k2*n)/(p_L))*(1-theta^2)*sin(2*u);

Drf_dot=((k2*n)/(p_L))*((1-theta^2)*cos(2*u)-(3/2)*(1-3*theta^2));

rk=r*(1-(3/2)*k2*(sqrt(1-e_L^2)/(p_L^2))*(3*theta^2-1))+Dr;

uk=u+Du;

Omegak=Omega+DOmega;
ik=io+Di;
r_dotk=r_dot+Dr_dot;
rf_dotk=rf_dot+Drf_dot;

M=[-sin(Omegak)*cos(ik);cos(Omegak)*cos(ik);sin(ik)];

N=[cos(Omegak);sin(Omegak);0];

U=M*sin(uk)+N*cos(uk);
V=M*cos(uk)-N*sin(uk);

r=radiusearthkm*rk*U;
v=vkmpersec*(r_dotk*U+(rf_dotk)*V);

oe_mean_out=[a;e;io;Omega;omega;u];

% update sat structure 

sat. ndot= sat. ndot;
sat.nddot= sat.nddot;

sat.inclo= io;
sat.nodeo= Omega;
sat.ecco= e;
sat.argpo= omega;
sat.mo= Mp;
sat.no= n;
sat.a= a;

sat.alta=a*(1+e)-1;
sat.altp=a*(1-e)-1;

sat.jdsatepoch=sat.jdsatepoch+(t-to)/(3600*24);
sat.error= 0;
sat.epochyr= 22;
sat.epochdays= 232.2671;




