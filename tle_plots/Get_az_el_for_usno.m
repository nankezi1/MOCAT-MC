close all 
clear all 
clc
% load('oe_data')

load('catolog_data_geos')
% save('catolog_data_geos','satreci','X_eci','X_ecf','longstr1','longstr2')
npts=length(X_eci(1,:,1));
numcat=length(X_eci(1,1,:));radiusearthkm =6378.137;
NOFS_long=-77.066945;NOFS_lat=38.919254;NOFS_alt=0.070/1000;
azm_opt=zeros(npts,numcat);elm_opt=zeros(npts,numcat);avail_opt=zeros(npts,numcat);
avail_opt_camera=zeros(npts,numcat,7);b_sight=zeros(3,npts,numcat);ym=zeros(2,npts,numcat,7);
for i=1:numcat,
 [az,el,~,x_sight]=az_el_range(X_eci(1:3,:,i)',NOFS_long,NOFS_lat,NOFS_alt,gst(:,:,i)'*180/pi);
 avail_opt(find(el > 5 & el < 90),i)=1;
 azm_opt(:,i)=az;elm_opt(:,i)=el;
 ymi(:,:,i)=radec(X_eci(1:3,:,i))';

end
load('/Users/richlinares/Desktop/Rich/grants/2014/usno-cots/code/toolboxes/eop-processing/currenteopfile/eop_data.mat')
addpath(genpath('/Users/richlinares/Desktop/Rich/grants/2014/usno-cots/code/toolboxes'))
jd_current=57134.131085+2400000.5+2/24;
[year, month, day, ~, ~, ~] = jd2date(jd_current)
tdays=jd_current-2451545;
jcent=tdays/36525;
gst_sec=67310.54841+(876600*3600+8640184.812866)*jcent+0.093104*jcent.^2-6.2e-6*jcent.^3;
gstj=rem(gst_sec,86400)/240;
i_neg=find(gstj<0);gstj(i_neg)=gstj(i_neg)+360;

%Earth oreintation parameters SCALLING
%  MJD Xpole Ypole UT1-UTC LOD Xsig Ysig UTsig LODsig Nr Nf Nt Xrt Yrt Xrtsig Yrtsig

%             (10**-6")       (0.1 usec)    (10**-6")     (0.1 usec)              (10**-6"/d)    (10**-6"/d)
EOP_Scales=[1,1e-6,1e-6,      1e-7,   1e-7,     1e-6,1e-6, 1e-7,  1e-7, 1, 1, 1,   1e-6,   1e-6,  1e-6, 1e-6];
% set ECI conversion parameters
indexi=find(iers_iau2000_data.jd_utc<jd_current);
index=indexi(end);
% Parameters determined for gps week
EOP=[iers_iau2000_data.pm_x_arcsec(index) iers_iau2000_data.pm_x_arcsec(index) iers_iau2000_data.pm_y_arcsec(index)...
    iers_iau2000_data.ut1_utc(index) iers_iau2000_data.lod_msec(index) iers_iau2000_data.pm_pred(index)...
    iers_iau2000_data.pm_pred(index)  iers_iau2000_data.ut1_utc_pred(index) 33 0 0 0 1253 842 27 42];
% set ECI conversion parameters
EOP = EOP_Scales.*EOP;lod=EOP(5);xp=EOP(2);yp=EOP(3);eqeterms=0;ddpsi=0;ddeps=0;

 RAj=178.716; DECj=-6.304;
 b=(35786+radiusearthkm)*[cos(RAj*pi/180)*cos(DECj*pi/180);
 sin(RAj*pi/180)*cos(DECj*pi/180);
 sin(DECj*pi/180)];
                            
%convert fron Jd to year month..etc.
    [year,month,day,hour,minu,sec,dayweek,dategreg] = julian2greg(jd_current);
    [ut1, tut1, jdut1, utc, tai, tt, ttt, jdtt, tdb, ttdb, jdtdb, tcg, jdtcg, tcb, jdtcb ] ...
        = convtime (year,month,day,hour,minu,sec, 0, 0, 0 );
[gps_p_eci_km,gps_v_eci_km,aeci,Aecef2eci] = ecef2eci(b,b,zeros(3,1),ttt,jd_current,lod,xp,yp,eqeterms,ddpsi,ddeps);
r=Aecef2eci'*b;
 [azjj,eljj,~,x_sightjj]=az_el_range(b',NOFS_long,NOFS_lat,NOFS_alt,gstj*180/pi);

eli=[];
azi=[];
I=find(oe_geos(:,3)*180/pi<5 & abs(35786+radiusearthkm-oe_geos(:,1))<100);temp=[];
load('GEO_NAMES')
figure(5)
for i=I'
 if avail_opt(:,i)'*avail_opt(:,i)~=0 & azm_opt(find(avail_opt(:,i)==1),i)>170 & azm_opt(find(avail_opt(:,i)==1),i)<187
     i
%  plot(azm_opt(find(avail_opt(:,i)==1),i),elm_opt(find(avail_opt(:,i)==1),i),'.b','MarkerSize',10); hold on;
 a2=azm_opt(find(avail_opt(:,i)==1),i);
 b2=elm_opt(find(avail_opt(:,i)==1),i) ;
 % display the results
  metric_string = sprintf('%2.2f',oe_geos(i,7));
  J=find(satnum==oe_geos(i,7));
  temp=[temp;oe_geos(i,7) i];
  % mark objects above the threshold with a black circle
  if ~isempty(find(avail_opt(:,i)==1))
     plot(a2(1)+.01,b2(1)+.01,'ko');
     text(a2(1)+.05,b2(1)+.05,satname{J},'Color','k',...
       'FontSize',10,'FontWeight','bold','Rotation',270);hold on;
   a2(1)
  end
grid on ;hold on;
xlabel('az (Deg)')
ylabel('el (Deg)')
% pause
eli=[eli elm_opt(find(avail_opt(:,i)==1),i)];
azi=[azi azm_opt(find(avail_opt(:,i)==1),i)];
 end
end
save('fov_sats_042215','temp')

eli=[];
azi=[];
I=find(oe_geos(:,3)*180/pi<75 & abs(35786+radiusearthkm-oe_geos(:,1))<100);temp=[];
figure(8)
for i=I'
 if avail_opt(:,i)'*avail_opt(:,i)~=0 & azm_opt(find(avail_opt(:,i)==1),i)>160 & azm_opt(find(avail_opt(:,i)==1),i)<195
     i
%  plot(azm_opt(find(avail_opt(:,i)==1),i),elm_opt(find(avail_opt(:,i)==1),i),'.b','MarkerSize',10); hold on;
 a2=azm_opt(find(avail_opt(:,i)==1),i);
 b2=elm_opt(find(avail_opt(:,i)==1),i) ;
 % display the results
  metric_string = sprintf('%2.2f',oe_geos(i,7));
  temp=[temp;oe_geos(i,7)];
  % mark objects above the threshold with a black circle
  if ~isempty(find(avail_opt(:,i)==1))
     plot(a2(1)+.01,b2(1)+.01,'ko');
     text(a2(1)+.05,b2(1)+.05,metric_string,'Color','k',...
       'FontSize',10,'FontWeight','bold','Rotation',90);hold on;
   a2(1)
  end
grid on ;hold on;
xlabel('az (Deg)')
ylabel('el (Deg)')
% pause
eli=[eli elm_opt(find(avail_opt(:,i)==1),i)];
azi=[azi azm_opt(find(avail_opt(:,i)==1),i)];
 end
end

eli=[];
azi=[];
I=find(oe_geos(:,3)*180/pi<.1 & abs(35786+radiusearthkm-oe_geos(:,1))<1);temp=[];
figure(5)
for i=I'
 if avail_opt(:,i)'*avail_opt(:,i)~=0 
     i
%  plot(azm_opt(find(avail_opt(:,i)==1),i),elm_opt(find(avail_opt(:,i)==1),i),'.b','MarkerSize',10); hold on;
 a2=azm_opt(find(avail_opt(:,i)==1),i);
 b2=elm_opt(find(avail_opt(:,i)==1),i) ;
 % display the results
  metric_string = sprintf('%2.2f',oe_geos(i,7));
  temp=[temp;oe_geos(i,7)];
  % mark objects above the threshold with a black circle
  if ~isempty(find(avail_opt(:,i)==1))
     plot(a2(1)+.01,b2(1)+.01,'ko');
%      text(a2(1)+.05,b2(1)+.05,metric_string,'Color','k',...
%        'FontSize',10,'FontWeight','bold','Rotation',90);hold on;
%    a2(1)
  end
grid on ;hold on;
xlabel('az (Deg)')
ylabel('el (Deg)')
% pause
eli=[eli elm_opt(find(avail_opt(:,i)==1),i)];
azi=[azi azm_opt(find(avail_opt(:,i)==1),i)];
 end
end



para_geo=polyfit(azi(find(eli > 5 & eli < 90)),eli(find(eli > 5 & eli < 90)),6);
y=polyval(para_geo,linspace(112,247,100));
% figure(190)
%  plot(linspace(112,247,100),y,'--b','MarkerSize',10)
% grid on ;hold on;
%  plot(azi(find(eli > 5 & eli < 90 & azi > 170 & azi <190)),eli(find(eli > 5 & eli < 90 & azi > 170 & azi <190)),'.r','MarkerSize',10)
% xlabel('az (Deg)')
% ylabel('el (Deg)')




azj=linspace(112+10,247-10,7);
azj(1)=azj(1)+5;azj(end)=azj(end)-5;
elj=polyval(para_geo,azj);
y=polyval(para_geo,linspace(112,247,100));
k = polyder(para_geo);dy=polyval(k,azj);
tilt=180/pi*tan(dy);
h=find(el > 5 & el < 90);
ii=floor(linspace(1,419,7));
cameras=cell(7,1);
for j=1:7
[line,az_fov,el_fov,o,azo(j),elo(j),R]=FOV_linesxy(10*pi/180,7.5*pi/180,azj(j)*pi/180,-elj(j)*pi/180,tilt(j)*pi/180);
camera.center_azel=[azo(j),elo(j)];
camera.center_azfov=az_fov;
camera.center_elfov=el_fov;
camera.line=line;
camera.center_unit_ecf=o;
camera.frame=R';
camera.theta=10*pi/180;
cameras{j}=camera;
 for i=1:4
     figure(200)
  plot(az_fov(:,i),el_fov(:,i),'.b','MarkerSize',10)
  hold on
  axis equal
  axis([110 245 5 60])
 end
 plot(azo(j),elo(j),'.r','MarkerSize',10)
end

     figure(200)
  plot(azjj,eljj,'*g','MarkerSize',20)


plot(linspace(112,247,100),y,'--m')
grid on ;hold on;
xlabel('az (Deg)')
ylabel('el (Deg)')


figure(100)
plot3(line(1,:),line(2,:),line(3,:),'.r')
hold on; grid on;
plot3(0,0,0)

obs=[];obsm=[];
for i=1:numcat,
 [az,el,~,x_sight]=az_el_range(X_eci(1:3,:,i)',NOFS_long,NOFS_lat,NOFS_alt,gst(:,:,i)'*180/pi);
 avail_opt(find(el > 5 & el < 90),i)=1;
 azm_opt(:,i)=az;elm_opt(:,i)=el;
 for k=1:npts
 b_sight(:,k,i)=x_sight(k,:)'/norm(x_sight(k,:));
 end
 for j=1:7
camera=cameras{j};
R=camera.frame;theta=camera.theta;b=R*b_sight(:,:,i);
% if b(1,:)<0.01
%     x=b(2,:);
%     y=b(3,:);
% else
x=b(2,:)./abs(b(1,:));y=b(3,:)./abs(b(1,:));
% end
avail_opt_camera(find(abs(x) < tan(theta) & abs(y) < tan(theta) & b(1,:)>0),i,j)=1;
if isempty(find(abs(x) < tan(theta) & abs(y) < tan(theta) & b(1,:)>0)) & ~isempty(find(el > 5 & el < 90))
   miss(i)=i; 
   obsm=[obsm; az(find(el > 5 & el < 90)) el(find(el > 5 & el < 90))];
end
ym(:,:,i,j)=radec(X_eci(1:3,:,i));

obs=[obs;az(find(avail_opt_camera(:,i,j)==1)) el(find(avail_opt_camera(:,i,j)==1))];

% figure(j)
% plot(x,y,'.r')
% hold on 
% axis([-tan(theta) tan(theta) -tan(theta) tan(theta)])
 end
end

figure
plot(obs(:,1),obs(:,2),'.r'); grid on;
xlabel('az (Deg)');ylabel('el (Deg)')

figure
plot(obsm(:,1),obsm(:,2),'.r'); grid on;hold on;
plot(obsj(:,1),obsj(:,2),'.b')
plot(obs(:,1),obs(:,2),'.g')
xlabel('az (Deg)');ylabel('el (Deg)')

