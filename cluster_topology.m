%���°�Ľ�Salam����������������㷨ͨ��MATLABԴ��
%{
������Ϊ���°�Դ�룬Դ����ɾ�����ܻ��Ư������������ͼƬ���㷨�Ľ�˵�����£�
1.ʹ��K��ֵ������ƽڵ�ֲ������ܣ�ʹ�ò���������������ͨ�Ժ;����Ը���
2.�����������������ݷḻ����������·�ķ��á�ʱ�ӡ������ڵ�ķ��á�ʱ�ӡ�ʱ�Ӷ�����������
3.��·ʱ�ӵ��ڽڵ�����������֮�����٣����ӷ���ʵ�����
%}

function [Sxy,AM,EdgeCost,EdgeDelay,EdgeBandWide,VertexCost,VertexDelay,VertexDelayJitter,VertexPacketLoss,ClusterMatrix]=cluster_topology(BorderLength,NodeAmount, ...
Alpha,Beta,PlotIf,EdgeCostDUB,EdgeBandWideDUB,VertexCostDUB,VertexDelayDUB,VertexDelayJitterDUB,VertexPacketLossDUB)
%% ��������б�
%BorderLength������������������ı߳�����λ��km
%NodeAmount������������ڵ�ĸ���
%Alpha����������������������AlphaԽ�󣬶̱���Գ��ߵı���Խ��
%Beta����������������������BetaԽ�󣬱ߵ��ܶ�Խ��
%PlotIf���������Ƿ���������ͼ�����Ϊ1����ͼ�����򲻻�ͼ
%EdgeCostDUB����������·���õĿ��Ʋ�����1*2���洢��·���õ��½���Ͻ�
%EdgeBandWideDUB����������·����Ŀ��Ʋ�����1*2���洢�½���Ͻ�
%VertexCostDUB���������ڵ���õĿ��Ʋ�����1*2,�洢�ڵ���õ��½���Ͻ�
%VertexDelayDUB���������ڵ�ʱ�ӵĿ��Ʋ�����1*2���ڴ��ڵ�ʱ�ӵ��½���Ͻ�
%VertexDelayJitterDUB���������ڵ�ʱ�Ӷ����Ŀ��Ʋ�����1*2���洢�ڵ�ʱ�Ӷ������½���Ͻ�
%VertexPacketLossDUB���������ڵ㶪���ʵĿ��Ʋ�����1*2,�洢�ڵ㶪���ʵ��½�

%%  �������
%Sxy��������3*N�ľ��󣬸��зֱ����ڴ洢�ڵ����ţ������꣬������ľ���
%AM��������0 1�洢����AM(i,j)=1��ʾ������i��j������ߣ�N*N
%EdgeCost����������·���þ���N*N
%EdgeDelay����������·ʱ�Ӿ���N*N
%EdgeBandWide����������·�������N*N 
%VertexCost���������ڵ��������,1*N
%VertexDelay���������ڵ�ʱ��������1*N
%VertexDelayJitter���������ڵ�ʱ�Ӷ�������,1*N
%VertexPacketLoss���������ڵ㶪��������,1*N
%%�Ƽ�������������� 
%BorderLength=1000;NodeAmount=25;Alpha=100000000;Beta=200000000000;
%PlotIf=1;EdgeCostDUB=[2,5];EdgeBandWideDUB=[30,1000];VertexCostDUB=[2,4];
%VertexDelayDUB=1e-4*[5,20];VertexDelayJitterDUB=1e-4*[3,8];
%VertexPacketLossDUB=1e-4*[0,500]
%%
t0 = clock;
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
% fprintf('RowCnt=%d ColCnt = %d RowBase = %f ColBase = %f\n',RowCnt,ColCnt,RowBase,ColBase);

%���������������������ѡȡNN���ڵ�
for i = 1:NN
    SSxy(i,1) = BorderLength*rand;
    SSxy(i,2) = BorderLength*rand;
end
%�ڵ���ȷֲ�
[IDX,C] = kmeans(SSxy,NodeAmount);Sxy = [[1:NodeAmount]',C]';


%�ִ�
for i = 1:NodeAmount
    row = floor(10 * Sxy(3,i)/(BorderLength*RowBase)) + 1;
    col = floor(10 * Sxy(2,i)/(BorderLength * ColBase)) + 1;
%     fprintf('row=%d col = %d x = %e y = %e\n',row,col,Sxy(2,i),Sxy(3,i));
    ClusterMatrix{row,col} = [ClusterMatrix{row,col},[Sxy(2,i),Sxy(3,i)]'];
    %���Դ��룬���ClusterMatrix
%     Cluster = ClusterMatrix{row,col};
%     [x,y] = size(Cluster);
%     fprintf('x=%d y = %d\n',x,y);
%     num = size(Cluster,2);
%     fprintf('num=%d\n',num);
%     for j = 1:num
%         fprintf('x=%e y = %e\n',Cluster(1,j),Cluster(2,j));
%     end
    
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
                error('������ϣ������ԣ�');
            end 
            tmp(1,k) = k;
            tmp(2,k) = Cluster(1,pos);
            tmp(3,k) = Cluster(2,pos);
        end
        ClusterMatrix{i,j} = tmp;
        %���Դ��룬���ClusterMatrix
%         Cluster = ClusterMatrix{i,j};
%         [x,y] = size(Cluster);
%         fprintf('x=%d y = %d\n',x,y);
%         num = size(Cluster,2);
%         fprintf('num=%d\n',num);
%         for m = 1:num
%             fprintf('idx = %d x=%e y = %e\n',Cluster(1,m),Cluster(2,m),Cluster(3,m));
%         end
    end
end


%���������ʼ������Ԫ���飩
cell(RowCnt,ColCnt);
AM = cell(RowCnt,ColCnt);
EdgeCost = cell(RowCnt,ColCnt);
EdgeDelay = cell(RowCnt,ColCnt);
EdgeBandWide = cell(RowCnt,ColCnt);
VertexCost  = cell(RowCnt,ColCnt);
VertexDelay = cell(RowCnt,ColCnt);
VertexDelayJitter = cell(RowCnt,ColCnt);
VertexPacketLoss  = cell(RowCnt,ColCnt);

for i = 1:RowCnt
    for j = 1:ColCnt
        Cluster = ClusterMatrix{i,j};
        nodes = size(Cluster,2);
        am = zeros(nodes,nodes);
        edgeCost = zeros(nodes,nodes);
        edgeDelay = zeros(nodes,nodes);
        edgeBandWide = zeros(nodes,nodes);
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
                    edgeBandWide(m,n) = EdgeBandWideDUB(1)+(EdgeBandWideDUB(2)-EdgeBandWideDUB(1))*rand;
                    edgeBandWide(n,m)=edgeBandWide(m,n);
                else
                    edgeDelay(m,n) = inf;
                    edgeDelay(n,m) = inf;
                    edgeCost(m,n) = inf;
                    edgeCost(n,m) = inf;
                    edgeBandWide(m,n) = inf;
                    edgeBandWide(n,m) = inf;
                end
            end
        end
        AM{i,j} = am;
        EdgeCost{i,j} = edgeCost;
        EdgeDelay{i,j} = edgeDelay;
        EdgeBandWide{i,j} = edgeBandWide;
    end
end

%����ͼ����ʾ��Χ
xlim([0,BorderLength]);
ylim([0,BorderLength]);
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
        VertexDelayJitter{i,j} = vertexDelayJitter;
        VertexPacketLoss{i,j} = vertexPacketLoss;
        Net_plot(BorderLength,nodes,Cluster,EdgeCost{i,j},PlotIf);
    end
end
if PlotIf == 1
    xlabel('x (km)','FontName','Times New Roman','FontSize',12);
    ylabel('y (km)','FontName','Times New Roman','FontSize',12);
end

timeDiff = etime(clock,t0);
fprintf('costTime = %f seconds\n',timeDiff);
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
if PlotIf == 1
    for i = 1:(nodes-1)
        for j = (i+1):nodes
            if isinf(edgeCost(i,j)) == 0
                plot([Cluster(2,i),Cluster(2,j)],[Cluster(3,i),Cluster(3,j)]);
                hold on;
            end
        end
    end
end
end