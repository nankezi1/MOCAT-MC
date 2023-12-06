close all 
clear all 
clc
% addpath(genpath('/Users/linares/Desktop/dod-dbase/grants/usno-cots/code/example code/valado'))
addpath(genpath('/Users/linares/Google Drive/research/Eye candy- general /tle_plots/vallado'))
addpath(genpath(pwd))

dt=1;npts=1441;tsince=0+[0:npts-1]*dt;
% filename='tle_single_27386.txt';
filename='tle_single_29106.txt';%Debris 
%filename='tle_single_00016.txt';
fid = fopen(filename);
line=1;
while feof(fid) == 0
    temp = fgetl(fid);
    line=line+1;
   %do something with line
end
fclose(fid)

numcat=floor(line/2);
global tumin radiusearthkm xke j2 j3 j4 j3oj2

% X_eci=zeros(6,npts,numcat);
% X_ecf=zeros(6,npts,numcat);
X_eci=zeros(6,numcat);
X_ecf=zeros(6,numcat);

oe_geos=zeros(numcat,8);
time=zeros(numcat,3);
long =zeros(npts,1,numcat);
gst =zeros(1,npts,numcat);
lat =zeros(npts,1,numcat);
height =zeros(npts,1,numcat);
satreci = cell(numcat,1);

jd_current=57134.131085+2400000.5+2/24;
jd=jd_current;
npts=1;

fid = fopen(filename);
figure(1)
for i=1:numcat
    longstr1{i} = fgets(fid);
    longstr2{i} = fgets(fid);
    if longstr1{i}==-1
        break
    end
% longstr1 =[temp1{1,i} ' ' temp2{1,i} ' ' temp3{1,i} '   ' ...
%     temp4{1,i} ' ' temp5{1,i} '  ' temp6{1,i} '  ' temp7{1,i} ' ' temp8{1,i} '  ' temp9{1,i}];
% longstr2 =[temp1{2,i} ' ' temp2{2,i} ' ' temp3{2,i} ' ' ...
%     temp4{2,i} ' ' temp5{2,i} ' ' temp6{2,i} ' ' temp7{2,i} ' ' temp8{2,i} ' ' temp9{2,i}];
[satrec, startmfe, stopmfe, deltamin] = twoline2rvMOD(longstr1{i},longstr2{i});
oe_geos(i,:)=[satrec.a*radiusearthkm; satrec.ecco;satrec.inclo;
              satrec.nodeo;satrec.argpo;satrec.mo;satrec.satnum;satrec.jdsatepoch];
          time(i,:)=[satrec.epochdays;satrec.epochyr;satrec.bstar];
% for n=1:npts
%    [satreci{i}, X_eci(:,n,i), X_ecf(:,n,i),gst(:,n,i)]=spg4_ecf(satrec,tsince(n)); %tsince(n) tsince(n)
% end

 [~, X_eci(:,i), X_ecf(:,i),~]=spg4_ecf(satrec,30);
end
fclose(fid);


NOFS_long=-111.74058;NOFS_lat=35.18419;NOFS_alt=2293.0/1000;
num_tle=floor(line/2);
[V,I]=sort(oe_geos(1:num_tle,8));
date=(oe_geos(I,8)-min(oe_geos(I,8)))/365;

figure
plot(date,oe_geos(I,1),'.r','MarkerSize',10)
grid on 
xlabel('year since epoch (years)')
ylabel('a (km)')

figure
plot(date,oe_geos(I,1).*(1+oe_geos(I,2)),'.r','MarkerSize',10)
grid on 
xlabel('year since epoch (years)')
ylabel('rp (km)')

mu=398600;
figure
plot(date,(oe_geos(I,1).*(1-oe_geos(I,2).^2)*mu).^0.5,'.r','MarkerSize',10)
grid on 
xlabel('year since epoch (years)')
ylabel('h (km)')

figure
plot(date,-mu*(2*oe_geos(I,1)).^-1,'.r','MarkerSize',10)
grid on 
xlabel('year since epoch (years)')
ylabel('E (km)')

figure
plot(date,oe_geos(I,2),'.r','MarkerSize',10)
grid on 
xlabel('year since epoch (years)')
ylabel('e ()')

figure
plot(date,oe_geos(I,3)*180/pi,'.r','MarkerSize',10)
grid on 
xlabel('year since epoch (years)')
ylabel('i (Deg)')

% figure
% plot((oe_geos(I,8)-min(oe_geos(1:9708,8)))/365,'.r','MarkerSize',10)


figure(3)
plot3(oe_geos(:,1)/radiusearthkm,oe_geos(:,2),oe_geos(:,3)*180/pi,'.r','MarkerSize',10)
grid on 
xlabel('a (Re)')
ylabel('e ()')
zlabel('i (Degs)')



figure

plot3(X_eci(1,:),X_eci(2,:),X_eci(3,:),'.r','MarkerSize',10)
grid on 
axis equal 

