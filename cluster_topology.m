%最新版改进Salam网络拓扑随机生成算法通用MATLAB源码
%{
本程序为最新版源码，源码无删减，能绘出漂亮的网络拓扑图片，算法改进说明如下：
1.使用K均值聚类控制节点分布的疏密，使得产生的网络拓扑连通性和均匀性更好
2.产生的网络拓扑数据丰富，包括：链路的费用、时延、带宽，节点的费用、时延、时延抖动、丢包率
3.链路时延等于节点距离除以三分之二光速，更加符合实际情况
%}

function [Sxy,AM,EdgeCost,EdgeDelay,EdgeBW,VertexCost,VertexDelay,VertexPacketLoss,ClusterMatrix]=cluster_topology(BorderLength,NodeAmount, ...
Alpha,Beta,PlotIf,EdgeCostDUB,EdgeBWDUB,VertexCostDUB,VertexDelayDUB,VertexDelayJitterDUB,VertexPacketLossDUB,VertexSpeedDUB,VertexDirecDUB)
%% 输入参数列表
%BorderLength————正方形区域的边长，单位：km
%NodeAmount————网络节点的个数
%Alpha————网络特征参数，Alpha越大，短边相对长边的比例越大
%Beta————网络特征参数，Beta越大，边的密度越大
%PlotIf————是否画网络拓扑图，如果为1，则画图，否则不画图
%EdgeCostDUB————链路费用的控制参数，1*2，存储链路费用的下界和上界
%EdgeBWDUB————链路带宽的控制参数，1*2，存储下界和上界
%VertexCostDUB————节点费用的控制参数，1*2,存储节点费用的下界和上界
%VertexDelayDUB————节点时延的控制参数，1*2，节储节点时延的下界和上界
%VertexDelayJitterDUB————节点时延抖动的控制参数，1*2，存储节点时延抖动的下界和上界
%VertexPacketLossDUB————节点丢包率的控制参数，1*2,存储节点丢包率的下界

%%  输出参数
%Sxy————3*N的矩阵，各列分别用于存储节点的序号，横坐标，纵坐标的矩阵
%AM————0 1存储矩阵，AM(i,j)=1表示存在由i到j的有向边，N*N
%EdgeCost————链路费用矩阵，N*N
%EdgeDelay————链路时延矩阵，N*N
%EdgeBW————链路带宽矩阵，N*N 
%VertexCost————节点费用向量,1*N
%VertexDelay————节点时延向量，1*N
%VertexPacketLoss————节点丢包率向量,1*N
%%推荐的输入参数设置 
%BorderLength=1000;NodeAmount=25;Alpha=100000000;Beta=200000000000;
%PlotIf=1;EdgeCostDUB=[2,5];EdgeBWDUB=[30,1000];VertexCostDUB=[2,4];
%VertexDelayDUB=1e-4*[5,20];VertexDelayJitterDUB=1e-4*[3,8];
%VertexPacketLossDUB=1e-4*[0,500]
%%
%参数初始化
NN = 10*NodeAmount;
SSxy = zeros(NN,2);

%簇参数
K = 5;
ClusterNodes = 2*K;    %簇内节点个数
Clusters = NodeAmount / ClusterNodes;   %簇数量
RowCnt = floor(sqrt(Clusters));
ColCnt = floor(Clusters / RowCnt);
ClusterMatrix = cell(RowCnt,ColCnt);
RowBase = 10/RowCnt;
ColBase = 10/ColCnt;
VDTime = [10,30]; %速度和方向更新间隔时间 TODO实际设为[60,300]
vdTime = VDTime(1) + (VDTime(2) - VDTime(1)) * rand;
%绘制动图参数
pic_num = 1;


%在正方形区域内随机均匀选取NN个节点
for i = 1:NN
    SSxy(i,1) = BorderLength*rand;
    SSxy(i,2) = BorderLength*rand;
end
%节点均匀分布
[IDX,C] = kmeans(SSxy,NodeAmount);

%添加节点移动属性（方向和速度）初始化
VertexSpeed = VertexSpeedDUB(1) + round(VertexSpeedDUB(2)*rand(NodeAmount,1));
VertexDirec = VertexDirecDUB(1) + round(VertexDirecDUB(2)*rand(NodeAmount,1));
Sxy = [[1:NodeAmount]',C,VertexSpeed,VertexDirec]';

%设置图形显示范围
xlim([0,BorderLength]);
ylim([0,BorderLength]);

%坐标名
if PlotIf == 1
    xlabel('x (km)','FontName','Times New Roman','FontSize',12);
    ylabel('y (km)','FontName','Times New Roman','FontSize',12);
end

 %输出参数初始化（胞元数组）
% global AM, EdgeCost,EdgeDelay,EdgeBW,VertexCost,VertexDelay,VertexPacketLoss;
AM = cell(RowCnt,ColCnt);
EdgeCost = cell(RowCnt,ColCnt);
EdgeDelay = cell(RowCnt,ColCnt);
EdgeBW = cell(RowCnt,ColCnt);
VertexCost  = cell(RowCnt,ColCnt);
VertexDelay = cell(RowCnt,ColCnt);
VertexPacketLoss  = cell(RowCnt,ColCnt);

%初始时间
t_pos = clock;
t_vd = clock;

%分簇
for i = 1:NodeAmount
    row = max(ceil(10 * Sxy(3,i)/(BorderLength*RowBase)),1);
    col = max(ceil(10 * Sxy(2,i)/(BorderLength * ColBase)),1);
    ClusterMatrix{row,col} = [ClusterMatrix{row,col},[Sxy(2,i),Sxy(3,i),Sxy(4,i),Sxy(5,i)]'];   
end

%为簇内节点编号
for i = 1:RowCnt
    for j = 1:ColCnt
        Cluster = ClusterMatrix{i,j};
        xCoordinate = Cluster(1,:);
        xCoordinate_sort = sort(xCoordinate);
        [rows,cols] = size(Cluster);
        tmp = zeros(rows+1,cols);
        for k = 1:cols  %结点数
            pos = find(xCoordinate == xCoordinate_sort(k));
            if length(pos)>1
                pos = pos(1);
            end 
            tmp(1,k) = k;
            tmp(2,k) = Cluster(1,pos);
            tmp(3,k) = Cluster(2,pos);
            tmp(4,k) = Cluster(3,pos);
            tmp(5,k) = Cluster(4,pos);
        end
        ClusterMatrix{i,j} = tmp;
    end
end

%设置链路的初始参数值（连接矩阵，EdgeDelay,EdgeCost,EdgeBW)
for i = 1:RowCnt
    for j = 1:ColCnt
        Cluster = ClusterMatrix{i,j};
        nodes = size(Cluster,2);
        am = zeros(nodes,nodes);
        edgeCost = zeros(nodes,nodes);
        edgeDelay = zeros(nodes,nodes);
        edgeBW = zeros(nodes,nodes);
        for m = 1:(nodes-1)
            for n = (m+1):nodes
                Distance = ( (Cluster(2,m)-Cluster(2,n))^2 +...
                    (Cluster(3,m)-Cluster(3,n))^2)^0.5;
                P = Beta*exp(-Distance^5/(Alpha*BorderLength));
                if P>rand
                    am(m,n) = 1;
                    am(n,m) = 1;
                    edgeDelay(m,n) = 0.5*Distance/100000;
                    edgeDelay(n,m) = edgeDelay(m,n);
                    edgeCost(m,n) = EdgeCostDUB(1)+(EdgeCostDUB(2)-EdgeCostDUB(1))*rand;
                    edgeCost(n,m)=edgeCost(m,n);
                    edgeBW(m,n) = EdgeBWDUB(1)+(EdgeBWDUB(2)-EdgeBWDUB(1))*rand;
                    edgeBW(n,m)=edgeBW(m,n);
                else
                    edgeDelay(m,n) = inf;
                    edgeDelay(n,m) = inf;
                    edgeCost(m,n) = inf;
                    edgeCost(n,m) = inf;
                    edgeBW(m,n) = inf;
                    edgeBW(n,m) = inf;
                end
            end
        end
        AM{i,j} = am;
        EdgeCost{i,j} = edgeCost;
        EdgeDelay{i,j} = edgeDelay;
        EdgeBW{i,j} = edgeBW;
    end
end       

%设置节点的初始属性值（VertexCost,VertexDelay,VetexDelayJitter,VertexPacketLoss)
%并画拓扑图
for i = 1:RowCnt
    for j = 1:ColCnt
        Cluster = ClusterMatrix{i,j};
        nodes = size(Cluster,2);
        vertexCost  = zeros(1,nodes);
        vertexDelay = zeros(1,nodes);
        vertexDelayJitter = zeros(1,nodes);
        vertexPacketLoss  = zeros(1,nodes);
        for k = 1:nodes
            vertexCost(k) = VertexCostDUB(1)+(VertexCostDUB(2)-VertexCostDUB(1))*rand;
            vertexDelay(k) = VertexDelayDUB(1)+(VertexDelayDUB(2)-VertexDelayDUB(1))*rand;
            vertexDelayJitter(k) = VertexDelayJitterDUB(1)+(VertexDelayJitterDUB(2)-VertexDelayJitterDUB(1))*rand;
            vertexPacketLoss(k) = VertexPacketLossDUB(1)+(VertexPacketLossDUB(2)-VertexPacketLossDUB(1))*rand;
        end
        VertexCost{i,j} = vertexCost;
        VertexDelay{i,j} = vertexDelay;
        VertexPacketLoss{i,j} = vertexPacketLoss;

        Net_plot(BorderLength,nodes,Cluster,EdgeCost{i,j},PlotIf);
    end
end

%模拟节点动态移动
while true
    %节点移动
    t_posDiff = etime(clock,t_pos);
    for i = 1:RowCnt
        for j = 1:ColCnt
            Cluster = ClusterMatrix{i,j};
            nodes = size(Cluster,2);
            k = 1;
            while k <= nodes
                fprintf('k = %d nodes = %d\n',k,nodes);
                Cluster(2,k) = Cluster(2,k) + Cluster(4,k)*cos(Cluster(5,k))*t_posDiff;
                Cluster(3,k) = Cluster(3,k) + Cluster(4,k)*sin(Cluster(5,k))*t_posDiff;
                %节点坐标越界
                if Cluster(2,k) < 0
                    Cluster(2,k) = 0;
                elseif Cluster(2,k) > BorderLength
                     Cluster(2,k) = BorderLength;
                end
                if Cluster(3,k) < 0
                    Cluster(3,k) = 0;
                elseif Cluster(3,k) > BorderLength
                    Cluster(3,k) = BorderLength;
                end
                %节点移动至其它簇(节点加入与退出）
                row = max(ceil(10 * Cluster(3,k)/(BorderLength*RowBase)),1);
                col = max(ceil(10 * Cluster(2,k)/(BorderLength*ColBase)),1);
                if (row ~= i) || (col ~= j)
%                     fprintf('row = %d i = %d col = %d j = %d\n',row,i,col,j);
                    [ClusterMatrix,AM,EdgeCost,EdgeDelay,EdgeBW,VertexCost,...
                        VertexDelay,VertexPacketLoss] = ...
                        VertexOutAndIn(ClusterMatrix,[i,j],k,Cluster(:,k),AM,...
                    EdgeCost,EdgeDelay,EdgeBW,VertexCost,VertexDelay,...
                    VertexPacketLoss,BorderLength,RowBase,ColBase,...
                    EdgeCostDUB,EdgeBWDUB,VertexCostDUB,VertexDelayDUB,...
                    VertexPacketLossDUB);
                    Cluster = ClusterMatrix{i,j};
                    nodes = size(Cluster,2);
%                 fprintf('Outnodes = %d\n',size(EdgeCost{i,j},2));
                else
                    ClusterMatrix{i,j} = Cluster;
                    k = k + 1;
                end
            end       
        end
    end
    t_pos = clock;
    
    %随机变换节点速度和方向
    t_vdDiff = etime(clock,t_vd);
    if t_vdDiff > vdTime
        for i = 1:RowCnt
            for j = 1:ColCnt
                Cluster = ClusterMatrix{i,j};
                nodes = size(Cluster,2);
                VertexSpeed = VertexSpeedDUB(1) + round(VertexSpeedDUB(2)*rand(nodes,1));
                VertexDirec = VertexDirecDUB(1) + round(VertexDirecDUB(2)*rand(nodes,1));
                Cluster(4,:) = VertexSpeed;
                Cluster(5,:) = VertexDirec;
                ClusterMatrix{i,j} = Cluster;
            end
        end
        vdTime = VDTime(1) + (VDTime(2)-VDTime(1)) * rand;
    end
    
    %绘制动图
    drawnow;
    %set(gcf,'Units','centimeter','Position',[5 5 50 50]); %设置图片大小
    F=getframe(gcf);
    I=frame2im(F);
    [I,map]=rgb2ind(I,256);
    if pic_num == 1
        imwrite(I,map,'E:\simulation\topology.gif','gif', 'Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(I,map,'E:\simulation\topology.gif','gif','WriteMode','append','DelayTime',0.2);
    end
    pic_num = pic_num + 1;
    hold off;
    
    for i = 1:RowCnt
        for j = 1:ColCnt
            Cluster = ClusterMatrix{i,j};
            nodes = size(Cluster,2);
%             fprintf('ClusterMatrix nodes = %d\n',nodes);
%             fprintf('EdgeCost nodes = %d\n',size(EdgeCost{i,j},2));
            Net_plot(BorderLength,nodes,Cluster,EdgeCost{i,j},PlotIf);
        end
    end
end
    
end

%用于绘制网络拓扑的函数
function Net_plot(BorderLength,nodes,Cluster,edgeCost,PlotIf)
%画节点
if PlotIf == 1
    plot(Cluster(2,:),Cluster(3,:),'ko','MarkerEdgeColor','b','MarkerFaceColor','g','MarkerSize',5);    %'ko'：黑色圆圈；'MarkerEdgeColor'：标记的边框颜色；'MarkerFaceColor'：标记的颜色
    hold on;    %表示可在原图上修改
    %节点标序号
    for i = 1:nodes
        Str = int2str(i);
        text(Cluster(2,i)+BorderLength/100,Cluster(3,i)+BorderLength/100,Str,'FontName','Times New Roman','FontSize',12);
        hold on;
    end
end
%画边
    for i = 1:(nodes-1)
        for j = (i+1):nodes
            if isinf(edgeCost(i,j)) == 0
                plot([Cluster(2,i),Cluster(2,j)],[Cluster(3,i),Cluster(3,j)]);
                hold on;
            end
        end
    end
end