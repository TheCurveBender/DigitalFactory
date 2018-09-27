function polyplot(p0,p1,p2,p3,colorL,colorR)
t=linspace(0,1,2);
l1=t.*p0'+(1-t).*p1';
l2=t.*p1'+(1-t).*p2';
l3=t.*p2'+(1-t).*p3';
lw=1.5;
plot(l1(1,:),l1(2,:),'-','color',colorL,'lineWidth',0.6*lw);
hold on
plot(l2(1,:),l2(2,:),'-','color',[0.7,0.7,0.7],'lineWidth',0.5*lw);
plot(l3(1,:),l3(2,:),'-','color',colorR,'lineWidth',0.6*lw);
axis equal
end