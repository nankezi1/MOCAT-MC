function load_GEO_sats_CARMHF(catnum,night)
cameras=['A' 'B' 'C' 'D' 'E' 'F' 'G'];
obs=[];ym=[0.0662182180411231 -0.00132638607501781];

for i=1:length(cameras)
    [A1,A2,A3,A4,A5,A6,A7,A8,A9,A10]=...
        textread(['xid.f' cameras(i) 'M'],...
        '%s %s %s %f %f %f %f %s %s %s', -1);
    for mm=1:length(catnum)
        I=find(A5==catnum(mm));
        if ~isempty(I)
            J=find(str2num(cell2mat(A2(I)))==str2num(night));
            for k=1:length(I)
                filelocation= ['f' cameras(i) night 'M.fix'];
                A = exist(filelocation);
                if A~=0 && ~isempty(J) && ~isempty(I) && ~isempty(A4)
                    AA=textread(filelocation);
                    for j=1:length(J)
                        I2=find(AA(:,1)==A4(I(J(j))));
                        obs=[obs;[catnum(mm)*ones(length(AA(I2,3)),1) AA(I2,3)+2400000.5-2430000.0 (AA(I2,6)-ym(1))*pi/180 (AA(I2,7)-ym(2))*pi/180]];
                    end
                end
            end
        end
    end
end
obsa=obs;
[obsi,ai,~]=unique(obsa,'rows');
    obsi=obsi(1:100,:);
    dlmwrite('/Users/richlinares/Desktop/Rich/grants/2014/usno-cots/CAR-MHF/CARMHF_usno_v1/car-mhf-matlab/SampleObs/pseudo_GEODSS.gen',obsi,'delimiter','\t','precision',10)
    