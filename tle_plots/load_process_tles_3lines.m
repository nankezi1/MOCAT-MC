function [oe_geos,X_eci,X_ecf,numcat,indx]=load_process_tles_3lines(name,time)


[A1,A2,A3,A4,A5,A6,A7,A8,A9]=textread(name,'%s %s %s %s %s %s %s %s %s',-1);
%addpath(genpath('/Users/richlinares/Desktop/Rich/grants/2014/usno-cots/code/toolboxes/vallado'))
% addpath(genpath('/Users/linares/Desktop/files/usno-cots/code/toolboxes/vallado'))
lines=length(A1)/3;
% temp1=reshape(A1,2,lines);
% temp2=reshape(A2,2,lines);
% satnums=temp2(2,:)';
% temp3=reshape(A3,2,lines);
% temp4=reshape(A4,2,lines);
% temp5=reshape(A5,2,lines);
% temp6=reshape(A6,2,lines);
% temp7=reshape(A7,2,lines);
% temp8=reshape(A8,2,lines);
% temp9=reshape(A9,2,lines);
dt=1;
% npts=1441;
npts=time;
tsince=0+[0:npts-1]*dt;
numcat=length(A1)/3;
global tumin radiusearthkm xke j2 j3 j4 j3oj2

X_eci=zeros(6,npts,numcat);
X_ecf=zeros(6,npts,numcat);
oe_geos=zeros(numcat,8);


jd_current=57134.131085+2400000.5+2/24;
jd=jd_current;
npts=1;
indx=[];

fid = fopen(name);
figure(1)
for i=1:numcat
    name=fgets(fid);
    longstr1{i} = fgets(fid);
    longstr2{i} = fgets(fid);
    
    
    k = strfind(name,'STARLINK');
    if ~isempty(k)
       indx=[indx;i]; 
    end
    [satrec, startmfe, stopmfe, deltamin] = twoline2rvMOD(longstr1{i},longstr2{i});
    oe_geos(i,:)=[satrec.a*radiusearthkm; satrec.ecco;satrec.inclo;
              satrec.nodeo;satrec.argpo;satrec.mo;satrec.satnum;satrec.jdsatepoch];
    for n=1:npts
        [satreci{i}, X_eci(:,n,i), X_ecf(:,n,i),gst(:,n,i)]=spg4_ecf(satrec,tsince(n)); %tsince(n) tsince(n)
    end
    i
end
fclose(fid);