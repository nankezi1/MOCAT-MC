% close all; 
% clear all;
% clc;
jd0=satrec.jdsatepoch;
jd=jd0+tsince/(24*60);

indx=floor(1441/6);

split=[1:indx:1441];
new_objects=zeros(length(split),1);total_obs=zeros(length(split),1);
for l=1:length(split)-1
 int=1;clear data
for j=1:7
for k=1:numcat

Y_meas=ym(:,find(avail_opt_camera(split(l):split(l+1),k,j)==1),k,j);
if ~isempty(Y_meas)
Y_meas=[Y_meas;jd(find(avail_opt_camera(split(l):split(l+1),k,j)==1))];
data{int,j}=Y_meas;
label(int)=k;
int=int+1;
end

end
end
split_data{l}=data;
total_obs(l+1)=int;
new_objects(l+1)=int-total_obs(l);
end


for i=1:length(split_data)
  data= split_data{i};
  [m,~]=size(data);
  for j=1:7
  for k=1:m
  Y_meas=data{k,j};
  plot(Y_meas(1,:),Y_meas(2,:))
  end
  end
end


figure 
plot(tsince(split),new_objects)
xlabel('Time (min)')
ylabel('New tracks')
grid on 