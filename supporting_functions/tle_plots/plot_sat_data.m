function plot_sat_data(catnum,num)
addpath('/Users/richlinares/Desktop/Rich/grants/2014/usno-cots/code/ukf-smooth-real')
        load('0924GL15')
bb=[242.2;-2.89]; %bais (arc-sec)
para.bb=bb;
jd_utc_start=X(1,1);
jd_utc_end=X(end,1);
t=(X(:,1)-X(1,1))*3600*24;
para.jd0=jd_utc_start;
para.t=t;
para.dt_dyn=1;
para.Am=1e-5;
[iers_data,nutation_data]=para_poly_fit(jd_utc_start,jd_utc_end);
para.nutation_data     = nutation_data;
para.iers_data         = iers_data;
para.fx=@forcemodel;
para.prop=@rk4prop;
para.hx=@observation;

%sight coordinates
lat = 35.18419;  %  wgs-84 latitude (deg)
lon = -111.74058;  %  wgs-84 longitude (deg)
alt = 2293.0/1000;  %  wgs-84 altitude (km)
%...conver lla to xyz fixed
[ lf(1), lf(2), lf(3) ] = WGS84_LLA_to_XYZ(lat, lon, alt);
para.lf=lf;
paraj=para;
cameras=['A' 'B' 'C' 'D' 'E' 'F' 'G'];
    for i=1:length(cameras)
        [A1,A2,A3,A4,A5,A6,A7,A8,A9,A10]=...
            textread(['xid.f' cameras(i) 'M'],...
            '%s %s %s %f %f %f %f %s %s %s', -1);
       I=find(A5==catnum);
       for k=1:length(I)
           k
       filelocation= ['f' cameras(i) A2{I(k)} 'M.fix'];
        A = exist(filelocation);
        if A~=0
        AA=textread(filelocation);
        I2=find(AA(:,1)==A4(I(k)));
        if ~isempty(I2)
            clear y
            int=0;
        for j=1:length(I2)

        if (X(1,1)<AA(I2(j),3)+2400000.5 && X(end,1)>AA(I2(j),3)+2400000.5 )    
          vq = interp1(X(:,1),X(:,2:end),AA(I2(j),3)+2400000.5)';
          paraj.tc=((AA(I2(j),3)+2400000.5)-jd_utc_start)*86400;
          y(:,j)= feval(paraj.hx,vq/1000,paraj);   
          int=int+1;
        end
        end
        end
        if int~=0
        figure(100)
        plot(y(1,:)*180/pi,y(2,:)*180/pi,'.b');hold on; 
        
        figure(num+20)
        plot(3600*(AA(I2,6)-y(1,:)'*180/pi),3600*(AA(I2,7)-y(2,:)'*180/pi),'.r','MarkerSize',20);
        xlabel('RA');
        ylabel('DEC');title(A10{I(1)});grid on;hold on  
        end
%         
%         figure(num)
%         plot(AA(I2,8),AA(I2,9),'.r');
%         hold on
%         cov(AA(I2,8))
%         cov(AA(I2,9))
%         
%         figure(num+1)
%         plot(AA(I2,3),AA(I2,10),'.r');
%         hold on
%         
%         figure(num+2)
%         plot(AA(I2,6),AA(I2,7),'.r');
%         hold on
% figure(num)
% xlabel('azimuth (Degs)');
% ylabel('elevation (Degs)');
% grid on
% title(A10{I(1)})
% hold on
% 
% figure(num+1)
% xlabel('Time (MJD)');
% ylabel('mag');
% set(gca,'YDir','reverse')
% title(A10{I(1)})
% grid on
% hold on
% 
% figure(num+2)
% xlabel('RA');
% ylabel('DEC');
% title(A10{I(1)})
% grid on
% hold on


%         pause ;        close all
        end
        end
    end
