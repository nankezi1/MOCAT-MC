close all; 
clear all; 
clc; 
[A1,A2,A3,A4,A5,A6,A7,A8,A9]=textread('geo_3line.txt','%s %s %s %s %s %s %s %s %s',-1);
addpath(genpath('/Users/richlinares/Desktop/Rich/grants/2014/usno-cots/code/toolboxes/vallado'))
lines=length(A1)/3;
temp1=reshape(A1,3,lines);
temp2=reshape(A2,3,lines);satnums=temp2(2,:)';
temp3=reshape(A3,3,lines);
temp4=reshape(A4,3,lines);
temp5=reshape(A5,3,lines);
temp6=reshape(A6,3,lines);
temp7=reshape(A7,3,lines);
temp8=reshape(A8,3,lines);
temp9=reshape(A9,3,lines);
dt=1;npts=1441;tsince=0+[0:npts-1]*dt;
numcat=length(A1)/3
fid = fopen('geo_3line.txt');
figure(1)
for i=1:numcat
    longstr1{i} = fgets(fid);
    longstr2{i} = fgets(fid);
    longstr3{i} = fgets(fid);
% longstr1 =[temp1{1,i} ' ' temp2{1,i} ' ' temp3{1,i} '   ' ...
%     temp4{1,i} ' ' temp5{1,i} '  ' temp6{1,i} '  ' temp7{1,i} ' ' temp8{1,i} '  ' temp9{1,i}];
% longstr2 =[temp1{2,i} ' ' temp2{2,i} ' ' temp3{2,i} ' ' ...
%     temp4{2,i} ' ' temp5{2,i} ' ' temp6{2,i} ' ' temp7{2,i} ' ' temp8{2,i} ' ' temp9{2,i}];
[satrec, startmfe, stopmfe, deltamin] = twoline2rvMOD(longstr2{i},longstr3{i});
satnum(i,1)=[satrec.satnum];
satname{i,1}=longstr1{i}(2:end);
end

save('GEO_NAMES''satnum','satname')