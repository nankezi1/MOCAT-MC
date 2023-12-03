close all 
clear all 
clc
global tumin radiusearthkm xke j2 j3 j4 j3oj2

time=1441;
npts=time;
addpath(genpath(pwd))

%name='well-tracked.txt';
%name='tle.txt';
name='tle_28884_2019.txt';
name='Starlink_tle_3line.txt';
%[oe_geos,X_eci,X_ecf,numcat]=load_process_tles(name,time);
[oe_geos,X_eci,X_ecf,numcat,indx]=load_process_tles_3lines(name,time);

long =zeros(npts,1,numcat);
gst =zeros(1,npts,numcat);
lat =zeros(npts,1,numcat);
height =zeros(npts,1,numcat);
satreci = cell(numcat,1);

oe_geos=oe_geos(indx,:);
X_ecf=X_ecf(:,:,indx);
numcat=length(indx);
% for i=1:numcat
% [yr, mth, day, hr, min, sec] = jd2date(satreci{i}.jdsatepoch);  
% [long(:,:,i),lat(:,:,i),height(:,:,i)]=gc2gd(X_eci(1:3,:,i)',yr,mth,day,hr,min,sec,60,24*3600,1);
% end

for i=1:1441
    figure(1)
    plot3(reshape(X_ecf(1,i,:),numcat,1,1),reshape(X_ecf(2,i,:),numcat,1,1),reshape(X_ecf(3,i,:),numcat,1,1),'.b','MarkerSize',10)
    hold on
    plot3(0,0,0)
    grid on
    xlabel('x (km)')
    ylabel('y (km)')
    zlabel('z (km)')
    axis equal
    
    
%     figure(2)
%     % Plot Results
%     plot(reshape(long(1,1,:),numcat,1,1),reshape(lat(1,1,:),numcat,1,1)+90,'.r','MarkerSize',20)
%     clear topo
end

%normal 
dates=datetime(oe_geos(:,end),'convertfrom','juliandate');
dates=datenum(dates);
figure;set(gcf,'Color','w');
plot(dates,oe_geos(:,3)*180/pi,'-k','LineWidth',2)
grid on;xlabel('Time (year)','FontSize',18);ylabel('INC (Deg)','FontSize',18); hold on;
plot(dates,oe_geos(:,3)*180/pi,'+k','MarkerSize',5)
tstart = datetime(2006,1,1,1,0,0);
tend = datetime(2019,1,1,1,0,0);
xlim([datenum(tstart) datenum(tend)]);
 ylim([-0.01 0.8])
set(gca,'FontSize',18);
datetick('x','yyyy','keeplimits')


%change
dates=datetime(oe_geos(:,end),'convertfrom','juliandate');
dates=datenum(dates);
figure;set(gcf,'Color','w');
plot(dates,oe_geos(:,3)*180/pi,'-b','LineWidth',2)
grid on;xlabel('Time (year)','FontSize',18);ylabel('i (Deg)','FontSize',18); hold on;
plot(dates,oe_geos(:,3)*180/pi,'.r','MarkerSize',5)
tstart = datetime(2010,1,1,1,0,0);
tend = datetime(2015,1,1,1,0,0);
xlim([datenum(tstart) datenum(tend)]);
ylim([-0 0.8])
set(gca,'FontSize',18);
datetick('x','yyyy','keeplimits')


%change (a)
dates=datetime(oe_geos(:,end),'convertfrom','juliandate');
dates=datenum(dates);
figure;set(gcf,'Color','w');
plot(dates,oe_geos(:,3),'+b','LineWidth',2)
grid on;xlabel('Time (year)','FontSize',18);ylabel('i (Deg)','FontSize',18); hold on;
% plot(dates,oe_geos(:,3)*180/pi,'.r','MarkerSize',5)
tstart = datetime(2013,0,0,0,0,0);
tend = datetime(2014,1,1,1,0,0);
xlim([datenum(tstart) datenum(tend)]);
% ylim([-0 0.8])
set(gca,'FontSize',18);
datetick('x','yyyy','keeplimits')


figure;
hist(oe_geos(find(oe_geos(:,1)-radiusearthkm<450),1)-radiusearthkm,1000)

figure(3)
plot3(oe_geos(:,1)/radiusearthkm,oe_geos(:,2),oe_geos(:,3)*180/pi,'.r','MarkerSize',5)
grid on 
xlabel('a (Re)')
ylabel('e ()')
zlabel('i (Degs)')

figure(4)
plot3(oe_geos(:,1).*(1+oe_geos(:,2)),oe_geos(:,1).*(1-oe_geos(:,2)),oe_geos(:,3)*180/pi,'.r','MarkerSize',5)
grid on 
xlabel('a (Re)')
ylabel('e ()')
zlabel('i (Degs)')
set(gca,'XScale','log')
set(gca,'YScale','log')



figure(5)
loglog(oe_geos(:,1).*(1+oe_geos(:,2)),oe_geos(:,1).*(1-oe_geos(:,2)),'.b','MarkerSize',5)
grid on 
%axis([6.4*10^3 10^5 6.4*10^3 4.7*10^4])
xlabel('Perigee Altitude (km)')
ylabel('Apogee Altitude (km)')


figure(40)
plot(oe_geos(:,4)*180/pi,oe_geos(:,6)*180/pi,'.r','MarkerSize',5)
grid on 
ylabel('Mo (Degs)')
xlabel('\Omega (Degs)')


figure
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