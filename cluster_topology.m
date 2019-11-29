%最新版改进Salam网络拓扑随机生成算法通用MATLAB源码
%{
本程序为最新版源码，源码无删减，能绘出漂亮的网络拓扑图片，算法改进说明如下：
1.使用K均值聚类控制节点分布的疏密，使得产生的网络拓扑连通性和均匀性更好
2.产生的网络拓扑数据丰富，包括：链路的费用、时延、带宽，节点的费用、时延、时延抖动、丢包率
3.链路时延等于节点距离除以三分之二光速，更加符合实际情况
%}

function [Sxy,AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay,ClusterMatrix]=cluster_topology(BorderLength,NodeAmount, ...
Alpha,Beta,PlotIf,EdgeCostDUB,EdgeBWDUB,VertexDelayDUB,VertexSpeedDUB,VertexDirecDUB)
%% 输入参数列表
%BorderLength————正方形区域的边长，单位：km
%NodeAmount————网络节点的个数
%Alpha————网络特征参数，Alpha越大，短边相对长边的比例越大
%Beta————网络特征参数，Beta越大，边的密度越大
%PlotIf————是否画网络拓扑图，如果为1，则画图，否则不画图
%EdgeCostDUB————链路费用的控制参数，1*2，存储链路费用的下界和上界
%EdgeBWDUB————链路带宽的控制参数，1*2，存储下界和上界
%VertexDelayDUB————节点时延的控制参数，1*2，节储节点时延的下界和上界

%%  输出参数
%Sxy————3*N的矩阵，各列分别用于存储节点的序号，横坐标，纵坐标的矩阵
%AM————0 1存储矩阵，AM(i,j)=1表示存在由i到j的有向边，N*N
%EdgeCost————链路费用矩阵，N*N
%EdgeDelay————链路时延矩阵，N*N
%EdgeBW————链路带宽矩阵，N*N 
%VertexDelay————节点时延向量，1*N
%%推荐的输入参数设置 
%BorderLength=1000;NodeAmount=25;Alpha=100000000;Beta=200000000000;
%PlotIf=1;EdgeCostDUB=[2,5];EdgeBWDUB=[30,1000];
%VertexDelayDUB=1e-4*[5,20];
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
VDTime = [5,10]; %速度和方向更新间隔时间 TODO实际设为[60,300]
vdTime = VDTime(1) + (VDTime(2) - VDTime(1)) * rand;
%绘制动图参数
pic_num = 1;

%链路有效长度
MaxLinkDistance = sqrt( (1/3*RowBase)^2 + (1/3*ColBase)^2 )*BorderLength/10;


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
AM = cell(RowCnt,ColCnt);
EdgeCost = cell(RowCnt,ColCnt);
EdgeDelay = cell(RowCnt,ColCnt);
EdgeBW = cell(RowCnt,ColCnt);
VertexDelay = cell(RowCnt,ColCnt);

%初始时间
t_pos = clock;
t_vd = clock;

%分簇
for i = 1:NodeAmount
    row = min(max(ceil(10 * Sxy(3,i)/(BorderLength*RowBase)),1),RowCnt);
    col = min(max(ceil(10 * Sxy(2,i)/(BorderLength * ColBase)),1),ColCnt);
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
        nodesNum = size(Cluster,2);
        am = zeros(nodesNum,nodesNum);
        edgeCost = zeros(nodesNum,nodesNum);
        edgeDelay = zeros(nodesNum,nodesNum);
        edgeBW = zeros(nodesNum,nodesNum);
        for m = 1:(nodesNum-1)
            for n = (m+1):nodesNum
                Distance = ( (Cluster(2,m)-Cluster(2,n))^2 +...
                    (Cluster(3,m)-Cluster(3,n))^2)^0.5;
                P = Beta*exp(-Distance^5/(Alpha*BorderLength));
                if P>rand && Distance < MaxLinkDistance
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

%设置节点的初始属性值（VertexDelay)
%并画拓扑图
for i = 1:RowCnt
    for j = 1:ColCnt
        Cluster = ClusterMatrix{i,j};
        nodesNum = size(Cluster,2);
        vertexDelay = zeros(1,nodesNum);
        for k = 1:nodesNum
            vertexDelay(k) = VertexDelayDUB(1)+(VertexDelayDUB(2)-VertexDelayDUB(1))*rand;
        end
        VertexDelay{i,j} = vertexDelay;

        Net_plot(BorderLength,nodesNum,Cluster,AM{i,j},PlotIf);
    end
end

%模拟节点动态移动
while true
    %节点移动
    t_posDiff = etime(clock,t_pos);
    for i = 1:RowCnt
        for j = 1:ColCnt
            Cluster = ClusterMatrix{i,j};
            nodesNum = size(Cluster,2);
            k = 1;
            while k <= nodesNum
                Cluster(2,k) = Cluster(2,k) + Cluster(4,k)*cos(Cluster(5,k))*t_posDiff;
                Cluster(3,k) = Cluster(3,k) + Cluster(4,k)*sin(Cluster(5,k))*t_posDiff;
                %节点坐标越界
                if Cluster(2,k) < 0
                    Cluster(2,k) = BorderLength * 0.1;
                elseif Cluster(2,k) > BorderLength
                     Cluster(2,k) = BorderLength * 0.9;
                end
                if Cluster(3,k) < 0
                    Cluster(3,k) = BorderLength * 0.1;
                elseif Cluster(3,k) > BorderLength
                    Cluster(3,k) = BorderLength * 0.9;
                end
                ClusterMatrix{i,j} = Cluster;
                %移除无效链路
                [AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay] = RemoveIneffectiveLink...
                    (MaxLinkDistance,k,nodesNum,ClusterMatrix,[i,j],AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay);
                %若有节点脱离，尝试重连
                [AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay]...
                    = ReConnect(MaxLinkDistance,ClusterMatrix,[i,j],k,AM,EdgeCost,EdgeDelay,...
                        EdgeBW,VertexDelay,EdgeCostDUB,EdgeBWDUB,VertexDelayDUB);
                
                %节点移动至其它簇(节点加入与退出）
                row = min(max(ceil(10 * Cluster(3,k)/(BorderLength*RowBase)),1),RowCnt);
                col = min(max(ceil(10 * Cluster(2,k)/(BorderLength*ColBase)),1),ColCnt);
                if (row ~= i) || (col ~= j)
                    [ClusterMatrix,AM,EdgeCost,EdgeDelay,EdgeBW,...
                    VertexDelay] = VertexOutAndIn...
                    (MaxLinkDistance,ClusterMatrix,[i,j],k,Cluster(:,k),AM,...
                    EdgeCost,EdgeDelay,EdgeBW,VertexDelay,...
                    BorderLength,RowBase,ColBase,RowCnt,ColCnt,...
                    EdgeCostDUB,EdgeBWDUB,VertexDelayDUB);
                    Cluster = ClusterMatrix{i,j};
                    nodesNum = size(Cluster,2);
                else
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
                nodesNum = size(Cluster,2);
                VertexSpeed = VertexSpeedDUB(1) + round(VertexSpeedDUB(2)*rand(nodesNum,1));
                VertexDirec = VertexDirecDUB(1) + round(VertexDirecDUB(2)*rand(nodesNum,1));
                Cluster(4,:) = VertexSpeed;
                Cluster(5,:) = VertexDirec;
                ClusterMatrix{i,j} = Cluster;
            end
        end
        vdTime = VDTime(1) + (VDTime(2)-VDTime(1)) * rand;
    end
                
    
    %绘制动图
    drawnow;
%     set(gcf,'Units','centimeter','Position',[5 5 50 50]); %设置图片大小
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
            nodesNum = size(Cluster,2);
            Net_plot(BorderLength,nodesNum,Cluster,AM{i,j},PlotIf);
        end
    end
end
    
end

%用于绘制网络拓扑的函数
function Net_plot(BorderLength,nodesNum,Cluster,am,PlotIf)
%画节点
if PlotIf == 1
    plot(Cluster(2,:),Cluster(3,:),'ko','MarkerEdgeColor','b','MarkerFaceColor','g','MarkerSize',5);    %'ko'：黑色圆圈；'MarkerEdgeColor'：标记的边框颜色；'MarkerFaceColor'：标记的颜色
    hold on;    %表示可在原图上修改
    %节点标序号
    for i = 1:nodesNum
        Str = int2str(i);
        text(Cluster(2,i)+BorderLength/100,Cluster(3,i)+BorderLength/100,Str,'FontName','Times New Roman','FontSize',12);
        hold on;
    end
end
%画边
    for i = 1:(nodesNum-1)
        for j = (i+1):nodesNum
            if am(i,j) == 1
                plot([Cluster(2,i),Cluster(2,j)],[Cluster(3,i),Cluster(3,j)]);
                hold on;
            end
        end
    end
end

%移除无效链路
function [AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay]...
    = RemoveIneffectiveLink(MaxLinkDistance,nodeIdx,nodesNum,ClusterMatrix_,ClusterIdx,AM_,EdgeCost_,EdgeDelay_,EdgeBW_,VertexDelay_)

    Cluster = ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)};
    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    
    %调试参数
    once = 1;
  
    for i = 1:nodesNum
        if am(nodeIdx,i) == 1 
            distance = ( (Cluster(2,nodeIdx)-Cluster(2,i))^2 +...
                    (Cluster(3,nodeIdx)-Cluster(3,i))^2)^0.5;
            if distance > MaxLinkDistance
                %调试
                if once == 1
                    fprintf('================================RemoveIneffectiveLink==================\n');
                    fprintf('clusterIdx=[%d,%d]\n',ClusterIdx(1),ClusterIdx(2));
                    for m = 1:nodesNum
                        for n = 1:nodesNum
                            fprintf('am(%d,%d) = %d\t',m,n,am(m,n));
                        end
                        fprintf('\n');
                    end
                end
                once = 2;
                fprintf('RemoveLink (%d,%d)\n',nodeIdx,i);
                [AM_,EdgeCost_,EdgeDelay_,EdgeBW_,VertexDelay_]...
                    = LinkDisconnect([nodeIdx,i],ClusterIdx,AM_,EdgeCost_,EdgeDelay_,EdgeBW_,VertexDelay_);
            end
        end
    end
    
    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    if once == 2
        for m = 1:nodesNum
            for n = 1:nodesNum
                fprintf('am(%d,%d) = %d\t',m,n,am(m,n));
            end
            fprintf('\n');
        end
        fprintf('==============================================================\n');
    end
    AM = AM_;
    EdgeCost = EdgeCost_;
    EdgeDelay = EdgeDelay_;
    EdgeBW = EdgeBW_;
    VertexDelay = VertexDelay_;
end

%链路断开
function [AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay]...
    = LinkDisconnect(LinkIdx,ClusterIdx,AM_,EdgeCost_,EdgeDelay_,EdgeBW_,VertexDelay_)

    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    edgeCost = EdgeCost_{ClusterIdx(1),ClusterIdx(2)};
    edgeDelay = EdgeDelay_{ClusterIdx(1),ClusterIdx(2)};
    edgeBW = EdgeBW_{ClusterIdx(1),ClusterIdx(2)};
    vertexDelay = VertexDelay_{ClusterIdx(1),ClusterIdx(2)};
    
    %更新AM   
    am(LinkIdx(1),LinkIdx(2)) = 0;
    am(LinkIdx(2),LinkIdx(1)) = 0;
    AM_{ClusterIdx(1),ClusterIdx(2)} = am; 
    
    
    %更新Edge*
    edgeCost(LinkIdx(1),LinkIdx(2)) = inf;
    edgeCost(LinkIdx(2),LinkIdx(1)) = inf;
    EdgeCost_{ClusterIdx(1),ClusterIdx(2)} = edgeCost;
    
    edgeDelay(LinkIdx(1),LinkIdx(2)) = inf;
    edgeDelay(LinkIdx(2),LinkIdx(1)) = inf;
    EdgeDelay_{ClusterIdx(1),ClusterIdx(2)} = edgeDelay;
    
    edgeBW(LinkIdx(1),LinkIdx(2)) = inf;
    edgeBW(LinkIdx(2),LinkIdx(1)) = inf;
    EdgeBW_{ClusterIdx(1),ClusterIdx(2)} = edgeBW;
    
    %更新Vertex*    
    vertexDelay(LinkIdx(1)) = inf;
    VertexDelay_{ClusterIdx(1),ClusterIdx(2)} = vertexDelay;
    
    AM = AM_;
    EdgeCost = EdgeCost_;
    EdgeDelay = EdgeDelay_;
    EdgeBW = EdgeBW_;
    VertexDelay = VertexDelay_;
end