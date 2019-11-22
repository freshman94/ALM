%���°�Ľ�Salam����������������㷨ͨ��MATLABԴ��
%{
������Ϊ���°�Դ�룬Դ����ɾ�����ܻ��Ư������������ͼƬ���㷨�Ľ�˵�����£�
1.ʹ��K��ֵ������ƽڵ�ֲ������ܣ�ʹ�ò���������������ͨ�Ժ;����Ը���
2.�����������������ݷḻ����������·�ķ��á�ʱ�ӡ������ڵ�ķ��á�ʱ�ӡ�ʱ�Ӷ�����������
3.��·ʱ�ӵ��ڽڵ�����������֮�����٣����ӷ���ʵ�����
%}

function [Sxy,AM,EdgeCost,EdgeDelay,EdgeBW,VertexCost,VertexDelay,VertexPacketLoss,ClusterMatrix]=cluster_topology(BorderLength,NodeAmount, ...
Alpha,Beta,PlotIf,EdgeCostDUB,EdgeBWDUB,VertexCostDUB,VertexDelayDUB,VertexDelayJitterDUB,VertexPacketLossDUB,VertexSpeedDUB,VertexDirecDUB)
%% ��������б�
%BorderLength������������������ı߳�����λ��km
%NodeAmount������������ڵ�ĸ���
%Alpha����������������������AlphaԽ�󣬶̱���Գ��ߵı���Խ��
%Beta����������������������BetaԽ�󣬱ߵ��ܶ�Խ��
%PlotIf���������Ƿ���������ͼ�����Ϊ1����ͼ�����򲻻�ͼ
%EdgeCostDUB����������·���õĿ��Ʋ�����1*2���洢��·���õ��½���Ͻ�
%EdgeBWDUB����������·����Ŀ��Ʋ�����1*2���洢�½���Ͻ�
%VertexCostDUB���������ڵ���õĿ��Ʋ�����1*2,�洢�ڵ���õ��½���Ͻ�
%VertexDelayDUB���������ڵ�ʱ�ӵĿ��Ʋ�����1*2���ڴ��ڵ�ʱ�ӵ��½���Ͻ�
%VertexDelayJitterDUB���������ڵ�ʱ�Ӷ����Ŀ��Ʋ�����1*2���洢�ڵ�ʱ�Ӷ������½���Ͻ�
%VertexPacketLossDUB���������ڵ㶪���ʵĿ��Ʋ�����1*2,�洢�ڵ㶪���ʵ��½�

%%  �������
%Sxy��������3*N�ľ��󣬸��зֱ����ڴ洢�ڵ����ţ������꣬������ľ���
%AM��������0 1�洢����AM(i,j)=1��ʾ������i��j������ߣ�N*N
%EdgeCost����������·���þ���N*N
%EdgeDelay����������·ʱ�Ӿ���N*N
%EdgeBW����������·�������N*N 
%VertexCost���������ڵ��������,1*N
%VertexDelay���������ڵ�ʱ��������1*N
%VertexPacketLoss���������ڵ㶪��������,1*N
%%�Ƽ�������������� 
%BorderLength=1000;NodeAmount=25;Alpha=100000000;Beta=200000000000;
%PlotIf=1;EdgeCostDUB=[2,5];EdgeBWDUB=[30,1000];VertexCostDUB=[2,4];
%VertexDelayDUB=1e-4*[5,20];VertexDelayJitterDUB=1e-4*[3,8];
%VertexPacketLossDUB=1e-4*[0,500]
%%
%������ʼ��
NN = 10*NodeAmount;
SSxy = zeros(NN,2);

%�ز���
K = 5;
ClusterNodes = 2*K;    %���ڽڵ����
Clusters = NodeAmount / ClusterNodes;   %������
RowCnt = floor(sqrt(Clusters));
ColCnt = floor(Clusters / RowCnt);
ClusterMatrix = cell(RowCnt,ColCnt);
RowBase = 10/RowCnt;
ColBase = 10/ColCnt;
VDTime = [10,30]; %�ٶȺͷ�����¼��ʱ�� TODOʵ����Ϊ[60,300]
vdTime = VDTime(1) + (VDTime(2) - VDTime(1)) * rand;
%���ƶ�ͼ����
pic_num = 1;


%���������������������ѡȡNN���ڵ�
for i = 1:NN
    SSxy(i,1) = BorderLength*rand;
    SSxy(i,2) = BorderLength*rand;
end
%�ڵ���ȷֲ�
[IDX,C] = kmeans(SSxy,NodeAmount);

%��ӽڵ��ƶ����ԣ�������ٶȣ���ʼ��
VertexSpeed = VertexSpeedDUB(1) + round(VertexSpeedDUB(2)*rand(NodeAmount,1));
VertexDirec = VertexDirecDUB(1) + round(VertexDirecDUB(2)*rand(NodeAmount,1));
Sxy = [[1:NodeAmount]',C,VertexSpeed,VertexDirec]';

%����ͼ����ʾ��Χ
xlim([0,BorderLength]);
ylim([0,BorderLength]);

%������
if PlotIf == 1
    xlabel('x (km)','FontName','Times New Roman','FontSize',12);
    ylabel('y (km)','FontName','Times New Roman','FontSize',12);
end

 %���������ʼ������Ԫ���飩
% global AM, EdgeCost,EdgeDelay,EdgeBW,VertexCost,VertexDelay,VertexPacketLoss;
AM = cell(RowCnt,ColCnt);
EdgeCost = cell(RowCnt,ColCnt);
EdgeDelay = cell(RowCnt,ColCnt);
EdgeBW = cell(RowCnt,ColCnt);
VertexCost  = cell(RowCnt,ColCnt);
VertexDelay = cell(RowCnt,ColCnt);
VertexPacketLoss  = cell(RowCnt,ColCnt);

%��ʼʱ��
t_pos = clock;
t_vd = clock;

%�ִ�
for i = 1:NodeAmount
    row = max(ceil(10 * Sxy(3,i)/(BorderLength*RowBase)),1);
    col = max(ceil(10 * Sxy(2,i)/(BorderLength * ColBase)),1);
    ClusterMatrix{row,col} = [ClusterMatrix{row,col},[Sxy(2,i),Sxy(3,i),Sxy(4,i),Sxy(5,i)]'];   
end

%Ϊ���ڽڵ���
for i = 1:RowCnt
    for j = 1:ColCnt
        Cluster = ClusterMatrix{i,j};
        xCoordinate = Cluster(1,:);
        xCoordinate_sort = sort(xCoordinate);
        [rows,cols] = size(Cluster);
        tmp = zeros(rows+1,cols);
        for k = 1:cols  %�����
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

%������·�ĳ�ʼ����ֵ�����Ӿ���EdgeDelay,EdgeCost,EdgeBW)
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

%���ýڵ�ĳ�ʼ����ֵ��VertexCost,VertexDelay,VetexDelayJitter,VertexPacketLoss)
%��������ͼ
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

%ģ��ڵ㶯̬�ƶ�
while true
    %�ڵ��ƶ�
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
                %�ڵ�����Խ��
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
                %�ڵ��ƶ���������(�ڵ�������˳���
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
    
    %����任�ڵ��ٶȺͷ���
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
    
    %���ƶ�ͼ
    drawnow;
    %set(gcf,'Units','centimeter','Position',[5 5 50 50]); %����ͼƬ��С
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

%���ڻ����������˵ĺ���
function Net_plot(BorderLength,nodes,Cluster,edgeCost,PlotIf)
%���ڵ�
if PlotIf == 1
    plot(Cluster(2,:),Cluster(3,:),'ko','MarkerEdgeColor','b','MarkerFaceColor','g','MarkerSize',5);    %'ko'����ɫԲȦ��'MarkerEdgeColor'����ǵı߿���ɫ��'MarkerFaceColor'����ǵ���ɫ
    hold on;    %��ʾ����ԭͼ���޸�
    %�ڵ�����
    for i = 1:nodes
        Str = int2str(i);
        text(Cluster(2,i)+BorderLength/100,Cluster(3,i)+BorderLength/100,Str,'FontName','Times New Roman','FontSize',12);
        hold on;
    end
end
%����
    for i = 1:(nodes-1)
        for j = (i+1):nodes
            if isinf(edgeCost(i,j)) == 0
                plot([Cluster(2,i),Cluster(2,j)],[Cluster(3,i),Cluster(3,j)]);
                hold on;
            end
        end
    end
end