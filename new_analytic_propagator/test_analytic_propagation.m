close all
clear all
clc;

re = 6378.137;
global mu req j2

mu = 398600.4415;
req = 6378.14;
j2 = 0.00018623;

%test function for method
% Martinusi, Vladimir, Lamberto Dell Elce, and Gaëtan Kerschen.
% "Analytic propagation of near-circular satellite orbits
% in the atmosphere of an oblate planet." Celestial Mechanics
% and Dynamical Astronomy 123, no. 1 (2015): 85-103.

H_0=300;
a_0=400+re;
e_0=0.2431;
inc_0=51*pi/180;
omega_0=pi/2;
bigO_0=45*pi/180;
Mo_0=0;

Cd=2.2;
S_ref_m=0.001;
rho_0=1e-13;

param.C0=(1/2*Cd*S_ref_m*1000);% in km units, 
param.Bstar=param.C0*2.461e-5;% use as Bstar, note we multiply by 2.461e-5

param.rho_0=rho_0;


t=linspace(0,24*3600*1,1000);
X=zeros(length(t),6);

X(1,:)=[a_0;e_0;inc_0;bigO_0;omega_0;Mo_0];
Y=zeros(length(t),6);

Y(1,:)=mean2osc_m(X(1,:)');

F107=(120-70)/2*cos(2*pi*(t/11/(24*3600*365)))+2*40+30*rand(1,length(t));
Ap=2;



for k=1:length(t)-1
    %osc to mean elements
    y_mean=X(k,:)';

    %set parameters
    param.F107=F107(k);
    param.Ap=Ap;
    param.t=t(k+1);
    param.t_0=t(k);

    %propagation mean element analytically
    [out_oe]=analytic_propagation(y_mean,param);

    %mean to osc elements
    %[y_mean]=mean_osculating_map(out_oe,1);

    %save outpus
    X(k+1,:)=out_oe';

    %Y(k+1,:) = mean_osculating_map(X(k+1,:)',-1);%from junkin's book
    Y(k+1,:) = mean2osc_m(X(k+1,:)');%from Matlabfile exchange
end



out1=mean_osculating_map(X(1,:)',1)
out2=mean_osculating_map(out1,-1)



figure;
subplot(3,2,1);plot(t/(24*3600*365),X(:,1)-re,'-r');grid on;
xlabel('Time (Day)');ylabel('a-re')

subplot(3,2,2);
plot(t/(24*3600),X(:,2),'-r');grid on;
xlabel('Time (Day)');ylabel('e')

subplot(3,2,3);
plot(t/(24*3600),X(:,3)*180/pi,'-r');grid on;
xlabel('Time (Day)');ylabel('inc (Deg)')

subplot(3,2,4);
plot(t/(24*3600),X(:,4)*180/pi,'-r');grid on;
xlabel('Time (Day)');ylabel('Omega (Deg)')

subplot(3,2,5);
plot(t/(24*3600),X(:,5)*180/pi,'-r');grid on;
xlabel('Time (Day)');ylabel('omega (Deg)')

subplot(3,2,6);
plot(t/(24*3600),X(:,6)*180/pi,'-r');grid on;
xlabel('Time (Day)');ylabel('Mo (Deg)')


%mean

figure;
subplot(3,2,1);plot(t/(24*3600*365),Y(:,1)-re,'-g');grid on;
xlabel('Time (Day)');ylabel('a-re')

subplot(3,2,2);
plot(t/(24*3600),Y(:,2),'-g');grid on;
xlabel('Time (Day)');ylabel('e')

subplot(3,2,3);
plot(t/(24*3600),real(Y(:,3))*180/pi,'-g');grid on;
xlabel('Time (Day)');ylabel('inc (Deg)')

subplot(3,2,4);
plot(t/(24*3600),Y(:,4)*180/pi,'-g');grid on;
xlabel('Time (Day)');ylabel('Omega (Deg)')

subplot(3,2,5);
plot(t/(24*3600),Y(:,5)*180/pi,'-g');grid on;
xlabel('Time (Day)');ylabel('omega (Deg)')

subplot(3,2,6);
plot(t/(24*3600),Y(:,6)*180/pi,'-g');grid on;
xlabel('Time (Day)');ylabel('Mo (Deg)')


