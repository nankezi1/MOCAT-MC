function [line,az,el,o,azo,elo,R]=FOV_lines(theta,az0,el0,tilt)
n=1000;
P0(:,1)=[1,0,-tan(theta)]';
P0(:,2)=[1,tan(theta),0]';
P0(:,3)=[1,0,tan(theta)]';
P0(:,4)=[1,-tan(theta),0]';

N(:,:,1)=[zeros(1,n);ones(1,n).*linspace(-tan(theta),tan(theta),n);zeros(1,n)];
N(:,:,2)=[zeros(1,n);zeros(1,n);ones(1,n).*linspace(-tan(theta),tan(theta),n)];
N(:,:,3)=[zeros(1,n);-ones(1,n).*linspace(-tan(theta),tan(theta),n);zeros(1,n)];
N(:,:,4)=[zeros(1,n);zeros(1,n);-ones(1,n).*linspace(-tan(theta),tan(theta),n)];

Rz=[cos(az0) -sin(az0) 0;sin(az0) cos(az0) 0;0 0 1];
% Ry=[cos(el0) 0 sin(el0);0 1 0;-sin(el0) 0 cos(el0)];
e=Rz*[0;1;0];
Ry=rod_formula(e,el0);
o=Ry*Rz*[1;0;0];
R=rod_formula(o,tilt)*Ry*Rz;
line=zeros(3,n,4);
% R=eye(3);
for i=1:4
line(:,:,i)=kron(P0(:,i),ones(1,n))+N(:,:,i);
for j=1:n
line(:,j,i)=R*line(:,j,i)/norm(line(:,j,i));
end
% Azimuth and Elevation of sight
az(:,i)=rem(atan2(line(2,:,i),line(1,:,i))*180/pi+360,360);
el(:,i)=asin(line(3,:,i))*180/pi;
end
o=R*[1;0;0];
azo=rem(atan2(o(2),o(1))*180/pi+360,360);
elo=asin(o(3))*180/pi;