%���°�Ľ�Salam����������������㷨ͨ��MATLABԴ��
%{
������Ϊ���°�Դ�룬Դ����ɾ�����ܻ��Ư������������ͼƬ���㷨�Ľ�˵�����£�
1.ʹ��K��ֵ������ƽڵ�ֲ������ܣ�ʹ�ò���������������ͨ�Ժ;����Ը���
2.�����������������ݷḻ����������·�ķ��á�ʱ�ӡ������ڵ�ķ��á�ʱ�ӡ�ʱ�Ӷ�����������
3.��·ʱ�ӵ��ڽڵ�����������֮�����٣����ӷ���ʵ�����
%}

function [Sxy,AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay,ClusterMatrix]=cluster_topology(BorderLength,NodeAmount, ...
Alpha,Beta,PlotIf,EdgeCostDUB,EdgeBWDUB,VertexDegreeDUB,VertexSpeedDUB,VertexDirecDUB)
%% ��������б�
%BorderLength������������������ı߳�����λ��km
%NodeAmount������������ڵ�ĸ���
%Alpha����������������������AlphaԽ�󣬶̱���Գ��ߵı���Խ��
%Beta����������������������BetaԽ�󣬱ߵ��ܶ�Խ��
%PlotIf���������Ƿ���������ͼ�����Ϊ1����ͼ�����򲻻�ͼ
%EdgeCostDUB����������·���õĿ��Ʋ�����1*2���洢��·���õ��½���Ͻ�
%EdgeBWDUB����������·����Ŀ��Ʋ�����1*2���洢�½���Ͻ�

%%  �������
%Sxy��������3*N�ľ��󣬸��зֱ����ڴ洢�ڵ����ţ������꣬������ľ���
%AM��������0 1�洢����AM(i,j)=1��ʾ������i��j������ߣ�N*N
%EdgeCost����������·���þ���N*N
%EdgeDelay����������·ʱ�Ӿ���N*N
%EdgeBW����������·�������N*N 
%VertexDelay���������ڵ�ʱ��������1*N
%%�Ƽ�������������� 
%BorderLength=1000;NodeAmount=25;Alpha=100000000;Beta=200000000000;
%PlotIf=1;EdgeCostDUB=[2,5];EdgeBWDUB=[30,1000];
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
VDTime = [5,10]; %�ٶȺͷ�����¼��ʱ�� TODOʵ����Ϊ[60,300]
vdTime = VDTime(1) + (VDTime(2) - VDTime(1)) * rand;
%���ƶ�ͼ����
pic_num = 1;

%��·��Ч����
MaxLinkDistance = sqrt( (1/3*RowBase)^2 + (1/3*ColBase)^2 )*BorderLength/10;


%���������������������ѡȡNN���ڵ�
for i = 1:NN
    SSxy(i,1) = BorderLength*rand;
    SSxy(i,2) = BorderLength*rand;
end
%�ڵ���ȷֲ�
[IDX,C] = kmeans(SSxy,NodeAmount);

%��ӽڵ��ƶ����ԣ�������ٶȣ���ʼ��
VertexSpeed = VertexSpeedDUB(1) + VertexSpeedDUB(2)*rand(NodeAmount,1);
VertexDirec = VertexDirecDUB(1) + VertexDirecDUB(2)*rand(NodeAmount,1);
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
AM = cell(RowCnt,ColCnt);
EdgeCost = cell(RowCnt,ColCnt);
EdgeDelay = cell(RowCnt,ColCnt);
EdgeBW = cell(RowCnt,ColCnt);
VertexDelay = cell(RowCnt,ColCnt);
VertexMaxDegree = cell(RowCnt,ColCnt);
LDT = cell(RowCnt,ColCnt);
VertexStability = cell(RowCnt,ColCnt);
VertexPriority = cell(RowCnt,ColCnt);
IsClusterHead = cell(RowCnt,ColCnt);

%��ʼʱ��
t_pos = clock;
t_vd = clock;

%�ִ�(Cluster��������Ϣ�У���ţ��ڵ����꣬�ڵ��ٶȣ��ڵ㷽��)
for i = 1:NodeAmount
    row = min(max(ceil(10 * Sxy(3,i)/(BorderLength*RowBase)),1),RowCnt);
    col = min(max(ceil(10 * Sxy(2,i)/(BorderLength * ColBase)),1),ColCnt);
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

%���ýڵ�����ȣ�VertexMaxDegree)
for i = 1:RowCnt
    for j = 1:ColCnt
        Cluster = ClusterMatrix{i,j};
        nodesNum = size(Cluster,2);    
        vertexMaxDegree = VertexDegreeDUB(1)+round( (VertexDegreeDUB(2)-VertexDegreeDUB(1)) * rand(1,nodesNum));
        VertexMaxDegree{i,j} = vertexMaxDegree;
    end
end


%������·�ĳ�ʼ����ֵ�����Ӿ���EdgeDelay,EdgeCost,EdgeBW)
for i = 1:RowCnt
    for j = 1:ColCnt
        Cluster = ClusterMatrix{i,j};
        nodesNum = size(Cluster,2);
        vertexMaxDegree = VertexMaxDegree{i,j};
        am = zeros(nodesNum,nodesNum);
        edgeCost = Inf(nodesNum,nodesNum);
        edgeDelay = Inf(nodesNum,nodesNum);
        edgeBW = Inf(nodesNum,nodesNum);
        ldt = zeros(nodesNum,nodesNum);
        vertexStability = zeros(1,nodesNum);
        for m = 1:(nodesNum-1)
            for n = (m+1):nodesNum
                Distance = ( (Cluster(2,m)-Cluster(2,n))^2 +...
                    (Cluster(3,m)-Cluster(3,n))^2)^0.5;
                P = Beta*exp(-Distance^5/(Alpha*BorderLength));
                %�ڵ��Լ��
                degree = sum(am(m,:));
                maxDegree = vertexMaxDegree(m);
                linkNodeDegree = sum(am(n,:));
                linkNodeMaxDegree = vertexMaxDegree(n);
                if degree >= maxDegree
                    break;
                end
                if linkNodeDegree >= linkNodeMaxDegree
                    continue;
                end
                if P>rand && Distance < MaxLinkDistance && degree < maxDegree
                    am(m,n) = 1;
                    am(n,m) = 1;
                    edgeDelay(m,n) = 0.5*Distance/100000;
                    edgeDelay(n,m) = edgeDelay(m,n);
                    edgeCost(m,n) = EdgeCostDUB(1)+(EdgeCostDUB(2)-EdgeCostDUB(1))*rand;
                    edgeCost(n,m)=edgeCost(m,n);
                    edgeBW(m,n) = EdgeBWDUB(1)+(EdgeBWDUB(2)-EdgeBWDUB(1))*rand;
                    edgeBW(n,m)=edgeBW(m,n);
                    ldt(m,n) = GetLDT(MaxLinkDistance,Cluster(:,[m,n]));
                    ldt(n,m) = ldt(m,n);
                end
            end
            vertexStability(m) = sum(ldt(m,:));
        end
        vertexStability(nodesNum) = sum(ldt(nodesNum,:));
        
        AM{i,j} = am;
        EdgeCost{i,j} = edgeCost;
        EdgeDelay{i,j} = edgeDelay;
        EdgeBW{i,j} = edgeBW;
        LDT{i,j} = ldt;
        VertexStability{i,j} = vertexStability;
    end
end


%���ýڵ�ĳ�ʼ����ֵ��VertexDelay)
for i = 1:RowCnt
    for j = 1:ColCnt
        Cluster = ClusterMatrix{i,j};
        nodesNum = size(Cluster,2);
        vertexDelay = zeros(1,nodesNum);
        edgeDelay = EdgeDelay{i,j};  
        am = AM{i,j};
        for m = 1:nodesNum
            vertexDelay(m) = GetVertexDelay(m,am,edgeDelay);
        end
        VertexDelay{i,j} = vertexDelay;
    end
end

%��ȡ�ڵ������ֵ
%��������ͼ
for i = 1:RowCnt
    for j = 1:ColCnt
        am = AM{i,j};
        vertexStability = VertexStability{i,j};
        vertexDelay = VertexDelay{i,j};
        vertexMaxDegree = VertexMaxDegree{i,j};
        nodesNum = size(vertexStability,2);
        vertexPriority = zeros(1,nodesNum);
        for m = 1:nodesNum
            vertexPriority(m) = GetPriority(m,vertexStability,am,...
                        vertexMaxDegree,vertexDelay);     
        end
        VertexPriority{i,j} = vertexPriority;
        IsClusterHead = SetClusterHead(IsClusterHead,[i,j],vertexPriority);
        Net_plot(BorderLength,nodesNum,ClusterMatrix{i,j},AM{i,j},PlotIf,...
            IsClusterHead{i,j},vertexPriority);
    end
end

%ģ��ڵ㶯̬�ƶ�
while true
    %�ڵ��ƶ�
    t_posDiff = etime(clock,t_pos);
    for i = 1:RowCnt
        for j = 1:ColCnt
            Cluster = ClusterMatrix{i,j};
            nodesNum = size(Cluster,2);
            k = 1;                    
            while k <= nodesNum
                Cluster(2,k) = Cluster(2,k) + Cluster(4,k)*cos(Cluster(5,k))*t_posDiff;
                Cluster(3,k) = Cluster(3,k) + Cluster(4,k)*sin(Cluster(5,k))*t_posDiff;
                %�ڵ�����Խ��
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
                
                %�Ƴ���Ч��·
                [IsClusterHead,AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay,LDT,VertexStability,...
                    VertexPriority] = RemoveIneffectiveLink(IsClusterHead,...
                    MaxLinkDistance,k,ClusterMatrix,[i,j],AM,EdgeCost,...
                    EdgeDelay, EdgeBW,VertexDelay,VertexMaxDegree,LDT,VertexStability,...
                    VertexPriority);
                %���нڵ����룬��������
                [IsClusterHead,AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay,LDT,...
                    VertexStability,VertexPriority] ...
                    = ReConnect(IsClusterHead,MaxLinkDistance,ClusterMatrix,...
                    [i,j],k,AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay,...
                    VertexMaxDegree,LDT,VertexStability,VertexPriority,...
                    EdgeCostDUB,EdgeBWDUB);
                
                %�ڵ��ƶ���������(�ڵ�������˳���
                row = min(max(ceil(10 * Cluster(3,k)/(BorderLength*RowBase)),1),RowCnt);
                col = min(max(ceil(10 * Cluster(2,k)/(BorderLength*ColBase)),1),ColCnt);
                if (row ~= i) || (col ~= j)
                    [IsClusterHead,ClusterMatrix,AM,EdgeCost,EdgeDelay,EdgeBW,...
                        VertexDelay,VertexMaxDegree,LDT,VertexStability,VertexPriority] ...
                        = VertexOutAndIn(IsClusterHead,MaxLinkDistance,ClusterMatrix,...
                        [i,j],k,Cluster(:,k),AM, EdgeCost,EdgeDelay,EdgeBW,VertexDelay,...
                        VertexMaxDegree,LDT,VertexStability,VertexPriority,...
                        BorderLength,RowBase,ColBase,RowCnt,ColCnt,...
                        EdgeCostDUB,EdgeBWDUB,VertexDegreeDUB);
                    Cluster = ClusterMatrix{i,j};
                    nodesNum = size(Cluster,2);
                else
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
                
    
    %���ƶ�ͼ
    drawnow;
%     set(gcf,'Units','centimeter','Position',[5 5 50 50]); %����ͼƬ��С
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
            Net_plot(BorderLength,nodesNum,Cluster,AM{i,j},PlotIf,IsClusterHead{i,j},VertexPriority{i,j});
        end
    end
end
    
end

%���ڻ����������˵ĺ���
function Net_plot(BorderLength,nodesNum,Cluster,am,PlotIf,isClusterHead,vertexPriority)
%���ڵ�
if PlotIf == 1
    plot(Cluster(2,:),Cluster(3,:),'ko','MarkerEdgeColor','b','MarkerFaceColor','g','MarkerSize',5);    %'ko'����ɫԲȦ��'MarkerEdgeColor'����ǵı߿���ɫ��'MarkerFaceColor'����ǵ���ɫ
    hold on;    %��ʾ����ԭͼ���޸�
    %�ڵ�����
    for i = 1:nodesNum
        priorityStr = sprintf('%.3f',vertexPriority(i));
        Str = [int2str(i),'(',priorityStr,')'];
        if isClusterHead(i) == 1
            plot(Cluster(2,i),Cluster(3,i),'r.','MarkerSize',30);
%             Str = strcat(Str,'(head)');
        end
        text(Cluster(2,i)+BorderLength/100,Cluster(3,i)+BorderLength/100,Str,'FontName','Times New Roman','FontSize',12);
        hold on;
    end
end
%����
    for i = 1:(nodesNum-1)
        for j = (i+1):nodesNum
            if am(i,j) == 1
                plot([Cluster(2,i),Cluster(2,j)],[Cluster(3,i),Cluster(3,j)]);
                hold on;
            end
        end
    end
end

%�Ƴ���Ч��·
function [IsClusterHead,AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay,LDT,VertexStability,...
    VertexPriority] = RemoveIneffectiveLink(IsClusterHead_,MaxLinkDistance,nodeIdx,...
    ClusterMatrix_,ClusterIdx,AM_,EdgeCost_,EdgeDelay_,EdgeBW_,VertexDelay_,...
    VertexMaxDegree_,LDT_,VertexStability_,VertexPriority_)

    Cluster = ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)};
    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    
    %���Բ���
    once = 1;
    nodesNum = size(am,2);
    for i = 1:nodesNum
        if am(nodeIdx,i) == 1 
            distance = ( (Cluster(2,nodeIdx)-Cluster(2,i))^2 +...
                    (Cluster(3,nodeIdx)-Cluster(3,i))^2)^0.5;
            if distance > MaxLinkDistance
                %����
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
                [IsClusterHead_,AM_,EdgeCost_,EdgeDelay_,EdgeBW_,VertexDelay_,LDT_,...
                    VertexStability_,VertexPriority_] = LinkDisconnect(IsClusterHead_,...
                    [nodeIdx,i],ClusterIdx,AM_,EdgeCost_,EdgeDelay_,EdgeBW_,...
                    VertexDelay_,VertexMaxDegree_,LDT_,VertexStability_,VertexPriority_);
            end
        end
    end
    IsClusterHead = IsClusterHead_;
    AM = AM_;
    EdgeCost = EdgeCost_;
    EdgeDelay = EdgeDelay_;
    EdgeBW = EdgeBW_;
    VertexDelay = VertexDelay_;
    LDT = LDT_;
    VertexStability = VertexStability_;
    VertexPriority = VertexPriority_;
end

%��·�Ͽ�
function [IsClusterHead,AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay,LDT,VertexStability,...
    VertexPriority] = LinkDisconnect(IsClusterHead_,LinkIdx,ClusterIdx,AM_,EdgeCost_,...
    EdgeDelay_,EdgeBW_,VertexDelay_,VertexMaxDegree_,LDT_,VertexStability_,...
    VertexPriority_)

    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    edgeCost = EdgeCost_{ClusterIdx(1),ClusterIdx(2)};
    edgeDelay = EdgeDelay_{ClusterIdx(1),ClusterIdx(2)};
    edgeBW = EdgeBW_{ClusterIdx(1),ClusterIdx(2)};
    vertexDelay = VertexDelay_{ClusterIdx(1),ClusterIdx(2)};
    ldt = LDT_{ClusterIdx(1),ClusterIdx(2)};
    vertexStability = VertexStability_{ClusterIdx(1),ClusterIdx(2)};
    vertexPriority = VertexPriority_{ClusterIdx(1),ClusterIdx(2)};
    vertexMaxDegree = VertexMaxDegree_{ClusterIdx(1),ClusterIdx(2)};
    
    %����AM   
    am(LinkIdx(1),LinkIdx(2)) = 0;
    am(LinkIdx(2),LinkIdx(1)) = 0;
    AM_{ClusterIdx(1),ClusterIdx(2)} = am; 
    
    %����Edge*
    edgeCost(LinkIdx(1),LinkIdx(2)) = inf;
    edgeCost(LinkIdx(2),LinkIdx(1)) = inf;
    EdgeCost_{ClusterIdx(1),ClusterIdx(2)} = edgeCost;
    
    edgeBW(LinkIdx(1),LinkIdx(2)) = inf;
    edgeBW(LinkIdx(2),LinkIdx(1)) = inf;
    EdgeBW_{ClusterIdx(1),ClusterIdx(2)} = edgeBW;
    
    edgeDelay(LinkIdx(1),LinkIdx(2)) = inf;
    edgeDelay(LinkIdx(2),LinkIdx(1)) = inf;
    EdgeDelay_{ClusterIdx(1),ClusterIdx(2)} = edgeDelay;
    
    %����LDT
    ldt(LinkIdx(1),LinkIdx(2)) = 0;
    ldt(LinkIdx(2),LinkIdx(1)) = 0;
    LDT_{ClusterIdx(1),ClusterIdx(2)} = ldt;
    
    %����Vertex*
    nodesNum = size(am,2);
    for i = 1:nodesNum
        vertexDelay(i) = GetVertexDelay(i,am,edgeDelay);
    end
    VertexDelay_{ClusterIdx(1),ClusterIdx(2)} = vertexDelay;
    
    for i = 1:nodesNum
        vertexPriority(i) = GetPriority(i,vertexStability,am,...
                vertexMaxDegree,vertexDelay);
    end
    VertexPriority_{ClusterIdx(1),ClusterIdx(2)} = vertexPriority;
    
    IsClusterHead = SetClusterHead(IsClusterHead_,ClusterIdx,vertexPriority);
    
    AM = AM_;
    EdgeCost = EdgeCost_;
    EdgeDelay = EdgeDelay_;
    EdgeBW = EdgeBW_;
    VertexDelay = VertexDelay_;
    LDT = LDT_;
    VertexStability = VertexStability_;
    VertexPriority = VertexPriority_;
end