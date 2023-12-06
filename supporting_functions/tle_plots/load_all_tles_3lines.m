function [sats,numcat,indx]=load_all_tles_3lines(name)


[A1,A2,A3,A4,A5,A6,A7,A8,A9] = textread(name,'%s %s %s %s %s %s %s %s %s',-1);

numcat=length(A1)/3;

X_eci=zeros(6,1,numcat);
X_ecf=zeros(6,1,numcat);

npts=1;
indx=[];

fid = fopen(name);
figure(1)

satreci=cell(numcat,1);
for i=1:numcat
    name=fgets(fid);
    longstr1{i} = fgets(fid);
    longstr2{i} = fgets(fid);
    [satrec, ~, ~, ~] = twoline2rvMOD(longstr1{i},longstr2{i});
    for n=1:npts
        [satreci{i}, X_eci(:,n,i), X_ecf(:,n,i),gst(:,n,i)]=spg4_ecf(satrec,0); %tsince(n) 
    end
    
end

sats=satreci;
fclose(fid);
