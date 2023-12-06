function y=radec(x)
r=(x(1,:).^2+x(2,:).^2+x(3,:).^2).^(0.5);
y=zeros(2,length(x(1,:)));
y(1,:)=atan2(x(2,:),x(1,:));
y(2,:)=asin(x(3)./r);