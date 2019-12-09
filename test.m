
clc;
clf;
[ClusterMatrix]=...
    cluster_topology(1e2,40,1e5,5e1,[3,5],[0,3],[0,2*pi]);
hold off;


% linkNodePos = [[1,5,5]',[2,10,10]'];
% speed = [1,2];
% direction = [1/6,1/6]*pi;
% MaxLinkDistance = 10;
% b = linkNodePos(2,1) - linkNodePos(2,2);
% d = linkNodePos(3,1) - linkNodePos(3,2);
% a = speed(1)*cos(direction(1)) - speed(2)*cos(direction(2));
% c = speed(1)*sin(direction(1)) - speed(2)*sin(direction(2));
%     
% LET = (-1*(a*b+c*d)+sqrt( (a.^2+ c.^2)*MaxLinkDistance.^2 - (a*d-b*c).^2))...
%             / (a.^2+ c.^2);
% fprintf('a = %f, b = %f, c = %f, d = %f, LET = %f\n',a,b,c,d,LET);

% pic_num = 1;
% for epsilon = 0.01:-0.001:0.005
%     t = 1;
%     syms x;
%     ur = -1;
%     ul = 1;
%     s = (ur + ul)/2;
%     w = ur + 1/2*(ul - ur)*(1-tanh((ul-ur)*(x-s*t)/(4*epsilon)));
%     figure(1);
%     ezplot(w);
%     axis([-0.05,0.05 -1.5 1.5])
%     drawnow;
%     F=getframe(gcf);
%     I=frame2im(F);
%     [I,map]=rgb2ind(I,256);
%     if pic_num == 1
%         imwrite(I,map,'E:\simulation\test.gif','gif', 'Loopcount',inf,'DelayTime',0.2);
%     else
%         imwrite(I,map,'E:\simulation\test.gif','gif','WriteMode','append','DelayTime',0.2);
%     end
%     pic_num = pic_num + 1;
% end
