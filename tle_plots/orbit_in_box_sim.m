close all; 
clear all; 
clc;

addpath(genpath(pwd))
load('catolog_data_geos')
tf=60*24*365*2;
dtj=1*60; % measurement sample time step
t=(0:dtj:tf); % time vector
nt=length(t);


i=520;
[satrec, startmfe, stopmfe, deltamin] = twoline2rvMOD(longstr1{i},longstr2{i});
NOFS_long=-111.74058;NOFS_lat=35.18419;NOFS_alt=2293.0/1000;


 [yr, mth, day, hr, min, sec] = jd2date(satrec.jdsatepoch);  

x_eci=zeros(6,nt);
gst=zeros(nt,1);
for i=1:nt
[satrecj, x_eci(:,i), ~,gst(i,1)]=spg4_ecf(satrec,t(i));
end


%  [az,el,alt,x_sight]=az_el_range(xx(1:3,:)',NOFS_long,NOFS_lat,NOFS_alt,gst*180/pi);
 [long,lat,height]=gc2gd(x_eci(1:3,:)',yr,mth,day,hr,min,sec,dtj*60,tf*60,1);
 
 figure 
 plot(long,lat,'.r')
 grid on 
 xlabel('Long (Deg)')
 ylabel('Lat (Deg)')
 
 
 figure
 plot3(long,lat,height)
 grid on 
 xlabel('Long (Deg)')
 ylabel('Lat (Deg)')
 zlabel('Height (km)')
 
 
 


figure(2)
load topo
imageHndl=imagesc(topo);
set(gca,'Ydir','normal','XLim',[0 360],'Ylim',[0 180])
colormap((topomap1+white)/2);
hold on;
set(gca,'Xtick',[0:30:360]')
set(gca,'XtickLabel',[0 30 60 90 120 150 180 -150 -120 -90 -60 -30 0])
set(gca,'Ytick',[0:10:180]')
set(gca,'YtickLabel',[-90,-80,-70,-60,-50,-40,-30,-20,-10,0, ...
10 20 30 40 50 60 70 80 90])
grid
xlabel('Longitude (Deg)')
ylabel('Latitude (Deg)')

plot(long,lat+90,'.r','MarkerSize',20)
