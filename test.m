%计算簇首脱离频率
% SpeedMax = 12;
% cnt = 5;
% HDFToS = zeros(cnt+1,SpeedMax);
% for i = 1:SpeedMax
%     HDFs = zeros(cnt,1);
%     for j = 1:cnt
%         [HDF]=cluster_topology(1e2,40,1e5,5e1,[3,5],[0,i],[0,2*pi]); 
%         HDFs(j) = HDF;
%     end
%     [HDFs;i]
%     HDFToS(:,i) = [HDFs;i];
% end
% HDFToS

%绘制簇首脱离频率图
% clc;
% clf;
% HDFToS = [0,0.5,1.2,2.6,3.2,4.2,6.5,7.2,10.4,11.75,13.2,15.2];
% C_HDFToS = [0.6,0.8,2.2,3.8,6.4,12,12.8,14.5,15.8,17.6,21.6,25];
% SpeedMax = 12;
% Speed = [1:SpeedMax];
% figure(2);
% hold on;
% box on;
% grid on;
% axis([0 14 0 28]);
% set(gca,'XTick',[1:1:14]);
% set(gca,'YTick',[0:2:28]);
% plot(Speed,HDFToS,'-ks','MarkerFaceColor','k');
% plot(Speed,C_HDFToS,'-ms','MarkerFaceColor','m');
% legend('MDD','DD','Location','southeast');
% hold off;

%绘制平均端到端时延图
% clc;
% clf;
% AP2PDToN = [12,13,14,14,15,15,15,16,17,17,17];
% C_AP2PDToN = [12,12,13,14,14,14,15,15,16,17,17];
% Node = [5:15];
% figure(3);
% hold on;
% box on;
% grid on;
% axis([4 16 11 18]);
% set(gca,'XTick',[5:1:16]);
% set(gca,'YTick',[12:1:18]);
% plot(Node,AP2PDToN,'-ks','MarkerFaceColor','k');
% plot(Node,C_AP2PDToN,'-ms','MarkerFaceColor','m');
% legend('MDD','DD','Location','southeast');
% hold off;

%计算平均端到端时延
% clc;
% clf;
% NodeInitial = 5;
% NodeMax = 15;
% cnt = 5;
% AP2PDToN = zeros(1,NodeMax-NodeInitial+1);
% for i = NodeInitial:NodeMax
%     delay = 0;
%     for j = 1:cnt
%         [AP2PD]=cluster_topology(i,1e3,4*i,1e9,5e1,[3,5],[0,2],[0,2*pi]);
%         AP2PD
%         delay = delay + AP2PD;
%     end
%     delay = delay/cnt;
%     i
%     delay 
%     AP2PDToN(i-4) = delay;
% end

clc;
clf;
[AP2PD]=cluster_topology(10,1e2,40,1e5,5e1,[3,5],[2,5],[0,2*pi]);

    
% A = [3,1;2,4;7,2;5,3];
% A
% B = sortrows(A,2);
% B
% B'

% a = cell(3,1);
% a{1} = [1,1];
% a{2} = [2,2,2];
% a{3} = [3,3,3,3];
% 
% a(2) = [];
% a
