close all 
clear all 
clc
[A1,A2,A3,A4,A5,A6,A7,A8,A9]=textread('tle.txt','%s %s %s %s %s %s %s %s %s',-1);
%addpath(genpath('/Users/richlinares/Desktop/Rich/grants/2014/usno-cots/code/toolboxes/vallado'))
addpath(genpath('/Users/linares/Desktop/files/usno-cots/code/toolboxes/vallado'))
lines=length(A1)/2;
temp1=reshape(A1,2,lines);
temp2=reshape(A2,2,lines);satnums=temp2(2,:)';
temp3=reshape(A3,2,lines);
temp4=reshape(A4,2,lines);
temp5=reshape(A5,2,lines);
temp6=reshape(A6,2,lines);
temp7=reshape(A7,2,lines);
temp8=reshape(A8,2,lines);
temp9=reshape(A9,2,lines);
dt=1;npts=1441;tsince=0+[0:npts-1]*dt;
numcat=length(A1)/2;
global tumin radiusearthkm xke j2 j3 j4 j3oj2

X_eci=zeros(6,npts,numcat);
X_ecf=zeros(6,npts,numcat);
oe_geos=zeros(numcat,7);
long =zeros(npts,1,numcat);
gst =zeros(1,npts,numcat);
lat =zeros(npts,1,numcat);
height =zeros(npts,1,numcat);
satreci = cell(numcat,1);

jd_current=57134.131085+2400000.5+2/24;
jd=jd_current;
npts=1;

fid = fopen('tle.txt');
figure(1)
for i=1:numcat
    longstr1{i} = fgets(fid);
    longstr2{i} = fgets(fid);
% longstr1 =[temp1{1,i} ' ' temp2{1,i} ' ' temp3{1,i} '   ' ...
%     temp4{1,i} ' ' temp5{1,i} '  ' temp6{1,i} '  ' temp7{1,i} ' ' temp8{1,i} '  ' temp9{1,i}];
% longstr2 =[temp1{2,i} ' ' temp2{2,i} ' ' temp3{2,i} ' ' ...
%     temp4{2,i} ' ' temp5{2,i} ' ' temp6{2,i} ' ' temp7{2,i} ' ' temp8{2,i} ' ' temp9{2,i}];
[satrec, startmfe, stopmfe, deltamin] = twoline2rvMOD(longstr1{i},longstr2{i});
oe_geos(i,:)=[satrec.a*radiusearthkm; satrec.ecco;satrec.inclo;
              satrec.nodeo;satrec.argpo;satrec.mo;satrec.satnum];
for n=1:npts
   [satreci{i}, X_eci(:,n,i), X_ecf(:,n,i),gst(:,n,i)]=spg4_ecf(satrec,tsince(n)); %tsince(n) tsince(n)
end


end
fclose(fid);

% save('catolog_data_geos_new','satreci','X_eci','X_ecf','longstr1','longstr2','gst','oe_geos')





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

% for i=1:numcat
% [yr, mth, day, hr, min, sec] = jd2date(satreci{i}.jdsatepoch);  
% [long(:,:,i),lat(:,:,i),height(:,:,i)]=gc2gd(X_eci(1:3,:,i)',yr,mth,day,hr,min,sec,60,24*3600,1);
% end

for i=1:1
    figure(1)
plot3(reshape(X_ecf(1,1,:),numcat,1,1),reshape(X_ecf(2,1,:),numcat,1,1),reshape(X_ecf(3,1,:),numcat,1,1),'.b','MarkerSize',10)
hold on 
plot3(0,0,0)
grid on 
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')
axis equal


figure(2)
% Plot Results
plot(reshape(long(1,1,:),numcat,1,1),reshape(lat(1,1,:),numcat,1,1)+90,'.r','MarkerSize',20)
clear topo
end

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
axis([6.4*10^3 10^5 6.4*10^3 4.7*10^4])
xlabel('Perigee Altitude (km)')
ylabel('Apogee Altitude (km)')


%     inc = str2num(tline(9:16));                                 % Orbit Inclination (degrees)
%     raan = str2num(tline(18:25));                               % Right Ascension of Ascending Node (degrees)
%     ecc = str2num(strcat('0.',tline(27:33)));                   % Eccentricity
%     w = str2num(tline(35:42));                                  % Argument of Perigee (degrees)
%     M = str2num(tline(44:51));                                  % Mean Anomaly (degrees)
%     sma = ( mu/(str2num(tline(53:63))*2*pi/86400)^2 )^(1/3);    % semi major axis
%     rNo = str2num(tline(64:68));                                % Revolution Number at Epoch 

% fid = fopen('tle.txt', 'rb');
% L1 = fscanf(fid,'%24c%*s',1);
% L2 = fscanf(fid,'%d%6d%*c%5d%*3c%*2f%f%f%5d%*c%*d%5d%*c%*d%d%5d',[1,9]);
% L3 = fscanf(fid,'%d%6d%f%f%f%f%f%f%f',[1,8]);
% 
%     load(fullfile(path2data,TLEfiles(2).name));
%     longstr1=line1;
%     longstr2=line2;
%     fprintf('USING 2-Line Elements: \n')
%     fprintf('%s \n',longstr1)
%     fprintf('%s \n',longstr2)
%     load(fullfile(path2data,TLEfiles(1).name));  %