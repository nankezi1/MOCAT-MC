% close all 
clear all 
clc
[A1,A2,A3,A4,A5,A6,A7,A8,A9]=textread('tle.txt','%s %s %s %s %s %s %s %s %s',-1);

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


fid = fopen('tle.txt');

satnum=[33275;32951;27854;28238;39122;33207;29494;27426];
numcati=length(satnum);
X_eci=zeros(6,npts,numcati);
X_ecf=zeros(6,npts,numcati);
oe_geos=zeros(numcati,6);
long =zeros(npts,1,numcati);
gst =zeros(1,npts,numcati);
lat =zeros(npts,1,numcati);
height =zeros(npts,1,numcati);
satreci = cell(numcati,1);
figure(1)

for i=1:numcat
    longstr1{i} = fgets(fid);
    longstr2{i} = fgets(fid);
% longstr1 =[temp1{1,i} ' ' temp2{1,i} ' ' temp3{1,i} '   ' ...
%     temp4{1,i} ' ' temp5{1,i} '  ' temp6{1,i} '  ' temp7{1,i} ' ' temp8{1,i} '  ' temp9{1,i}];
% longstr2 =[temp1{2,i} ' ' temp2{2,i} ' ' temp3{2,i} ' ' ...
%     temp4{2,i} ' ' temp5{2,i} ' ' temp6{2,i} ' ' temp7{2,i} ' ' temp8{2,i} ' ' temp9{2,i}];


[satrec, startmfe, stopmfe, deltamin] = twoline2rvMOD(longstr1{i},longstr2{i});
I=find(satnum==satrec.satnum);
if ~isempty(I)
    oe_geos(I,:)=[satrec.a*radiusearthkm; satrec.ecco;satrec.inclo;
              satrec.nodeo;satrec.argpo;satrec.mo];
for n=1:npts
   [satreci{I}, X_eci(:,n,I), X_ecf(:,n,I),gst(:,n,I)]=spg4_ecf(satrec,tsince(n));
end

end

end
fclose(fid);

% save('catolog_data_geos','satreci','X_eci','X_ecf','longstr1','longstr2')
%NOFS_long=-77.066945;NOFS_lat=38.919254;NOFS_alt=0.070/1000;
 NOFS_long=-111.74058;NOFS_lat=35.18419;NOFS_alt=2293.0/1000;


% Get Measurements
NOFS_min_el=20;NOFS_max_el=90;
azm_opt=zeros(npts,numcat);elm_opt=zeros(npts,numcat);avail_opt=zeros(npts,numcat);
avail_opt_camera=zeros(npts,numcat,7);b_sight=zeros(3,npts,numcat);ym=zeros(2,npts,numcat,7);
for i=1:numcati,
 [az,el,~,x_sight]=az_el_range(X_eci(1:3,:,i)',NOFS_long,NOFS_lat,NOFS_alt,gst(:,:,i)'*180/pi);
 avail_opt(find(el > NOFS_min_el & el < NOFS_max_el),i)=1;
 i
 figure(100);
 plot(az(find(el > NOFS_min_el & el < NOFS_max_el))-180,el(find(el > NOFS_min_el & el < NOFS_max_el)),'.r','MarkerSize',10)
xlabel('az (Deg) (South=0)');pause; grid on;
ylabel('el (Deg)')


figure(200);
 plot(az-180,el,'.r','MarkerSize',10)
xlabel('az (Deg) (South=0)');hold on; grid on;
ylabel('el (Deg)')
 
end

