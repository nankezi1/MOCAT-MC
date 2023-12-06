close all 
clear all 
clc
filename='obcom.txt';
[A1,A2,A3,A4,A5,A6,A7,A8,A9]=textread(filename,'%s %s %s %s %s %s %s %s %s',-1);

lines=length(A1)/2;
% temp1=reshape(A1,2,lines);
% temp2=reshape(A2,2,lines);satnums=temp2(2,:)';
% temp3=reshape(A3,2,lines);
% temp4=reshape(A4,2,lines);
% temp5=reshape(A5,2,lines);
% temp6=reshape(A6,2,lines);
% temp7=reshape(A7,2,lines);
% temp8=reshape(A8,2,lines);
% temp9=reshape(A9,2,lines);
dt=1;npts=60*60;tsince=0+[0:npts-1]*dt;
numcat=length(A1)/2;
global tumin radiusearthkm xke j2 j3 j4 j3oj2


fid = fopen(filename);

satnum=[33275;32951;27854;28238;39122;33207;29494;27426];
numcati=numcat;%length(satnum);
X_eci=zeros(6,npts,numcati);
X_ecf=zeros(6,npts,numcati);
oe_geos=zeros(numcati,7);
long =zeros(npts,1,numcati);
gst =zeros(1,npts,numcati);
lat =zeros(npts,1,numcati);
height =zeros(npts,1,numcati);
satreci = cell(numcati,1);
figure(1)
J0=2457024.79260538;
for i=1:numcat
    longstr1{i} = fgets(fid);
    longstr2{i} = fgets(fid);
% longstr1 =[temp1{1,i} ' ' temp2{1,i} ' ' temp3{1,i} '   ' ...
%     temp4{1,i} ' ' temp5{1,i} '  ' temp6{1,i} '  ' temp7{1,i} ' ' temp8{1,i} '  ' temp9{1,i}];
% longstr2 =[temp1{2,i} ' ' temp2{2,i} ' ' temp3{2,i} ' ' ...
%     temp4{2,i} ' ' temp5{2,i} ' ' temp6{2,i} ' ' temp7{2,i} ' ' temp8{2,i} ' ' temp9{2,i}];


[satrec, startmfe, stopmfe, deltamin] = twoline2rvMOD(longstr1{i},longstr2{i});
I=find(satnum==satrec.satnum);
% if ~isempty(I)
    oe_geos(i,:)=[satrec.a*radiusearthkm; satrec.ecco;satrec.inclo;
              satrec.nodeo;satrec.argpo;satrec.mo;satrec.jdsatepoch];
for n=1:npts
   [satreci{i}, X_eci(:,n,i), X_ecf(:,n,i),gst(:,n,i)]=spg4_ecf(satrec,(satrec.jdsatepoch-J0)*24*60+tsince(n)/60);
   [~,~,in_shadow(n,i), rsun(:,n,i),~,~] = sun_position(J0+tsince(n)/(24*3600), X_eci(1:3,n,i));
   A=5*2/(1000^2);   
end
vis_sol_const = 455*1000^2;

% end

end
fclose(fid);

% save('catolog_data_geos','satreci','X_eci','X_ecf','longstr1','longstr2')
%NOFS_long=-77.066945;NOFS_lat=38.919254;NOFS_alt=0.070/1000;
 NOFS_long=-111.74058;NOFS_lat=35.18419;NOFS_alt=2293.0/1000;
 
% Greenwich Mean Sidereal Time from Vallado (third edition), p. 191
jdate=2457024.79260538+tsince/60/(24*3600);
tdays=jdate-2451545;
jcent=tdays/36525;
gst_sec=67310.54841+(876600*3600+8640184.812866)*jcent+0.093104*jcent.^2-6.2e-6*jcent.^3;
gst=rem(gst_sec,86400)/240;
i_neg=find(gst<0);gst(i_neg)=gst(i_neg)+360;


% Get Measurements
NOFS_min_el=5;NOFS_max_el=90;
azm_opt=zeros(npts,numcat);elm_opt=zeros(npts,numcat);avail_opt=zeros(npts,numcat);
avail_opt_camera=zeros(npts,numcat,7);b_sight=zeros(3,npts,numcat);ym=zeros(2,npts,numcat,7);
for i=1:numcati,
 [az(:,i),el(:,i),~,x_sight]=az_el_range(X_eci(1:3,:,i)',NOFS_long,NOFS_lat,NOFS_alt,gst');
for j=1:length(x_sight(:,1))
    xLOS =(x_sight(j,:)-X_eci(1:3,j,i)')'/norm(x_sight(j,:)-X_eci(1:3,j,i)');
SCrad = 2/1000; % assuming a 2m radius for the spacecraft
k = vis_sol_const/((norm(X_eci(1:3,:,i)'))^2); % not sure where this came from
A = pi*SCrad^2; % area of the object
Rdif = 0.3; % assume completely difuse object and all light is reflected difusely
Flux = k*Rdif*A*(1+xLOS'*rsun(:,j,i)/norm(rsun(:,n,i)));
if in_shadow(j,i)==0
    m(j,i)=20;
else
m(j,i) = -26.7-2.5*log10(abs(Flux/vis_sol_const));
end
end

%  avail_opt(find(el > NOFS_min_el & el < NOFS_max_el),i)=1;
%  i
%  figure(100);
%  plot(az(find(el > NOFS_min_el & el < NOFS_max_el))-180,el(find(el > NOFS_min_el & el < NOFS_max_el)),'.r','MarkerSize',10)
% xlabel('az (Deg) (South=0)');hold on; grid on;
% ylabel('el (Deg)')

% figure(1000);
% plot(tsince,m(:,i));hold on;
% xlabel('az (Deg) (South=0)');hold on; grid on;
% ylabel('Apparent Magnitude')
% figure(200);
%  plot(az-180,el,'.r','MarkerSize',10)
% xlabel('az (Deg) (South=0)');hold on; grid on;
% ylabel('el (Deg)')
 
end

figure;
for i=1:numcati
plot(tsince(find(el(:,i) > NOFS_min_el & el(:,i) < NOFS_max_el))/60,m(find(el(:,i) > NOFS_min_el & el(:,i) < NOFS_max_el),i),'--b');hold on;
end
xlabel('Time (Min)');hold on; grid on;
ylabel('Apparent Magnitude')

% figure;
% for j=1:3600
% for i=1:numcati
% plot3(X_eci(1,j,i),X_eci(2,j,i),X_eci(3,j,i),'.-r','MarkerSize',10)
% hold on;
% pause(0.25)
% end
% end
figure 
for i=1:3600
plot(reshape(X_ecf(1,i,:),1,numcati),'*r'); 
pause(.01)
end
figure;
for i=1:numcati
plot(az(find(el(:,i) > NOFS_min_el & el(:,i) < NOFS_max_el),i)-180,el(find(el(:,i) > NOFS_min_el & el(:,i) < NOFS_max_el),i),'.r','MarkerSize',10)
hold on;
end
xlabel('az (Deg) (South=0)');hold on; grid on;
ylabel('el (Deg)')


figure;
for i=1:3600
    vis(i)=sum(el(i,:) > NOFS_min_el & el(i,:) < NOFS_max_el & in_shadow(i,:)~=0);
% plot(az(i,find(el(i,:) > NOFS_min_el & el(i,:) < NOFS_max_el)),el(i,find(el(i,:) > NOFS_min_el & el(i,:) < NOFS_max_el)),'.r','MarkerSize',10)
% ;hold on;
end
plot(tsince/60,vis)
xlabel('Time (Min)');hold on; grid on;
ylabel('Number of Satellites Vissible')
axis([0 60 0.5 3.5])

figure;
for i=1:numcati
plot(in_shadow(find(el(:,i) > NOFS_min_el & el(:,i) < NOFS_max_el),i),'-r');hold on 
end