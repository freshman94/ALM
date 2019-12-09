%���°�Ľ�Salam����������������㷨ͨ��MATLABԴ��
%{
������Ϊ���°�Դ�룬Դ����ɾ�����ܻ��Ư������������ͼƬ���㷨�Ľ�˵�����£�
1.ʹ��K��ֵ������ƽڵ�ֲ������ܣ�ʹ�ò���������������ͨ�Ժ;����Ը���
2.�����������������ݷḻ����������·�ķ��á�ʱ�ӡ������ڵ�ķ��á�ʱ�ӡ�ʱ�Ӷ�����������
3.��·ʱ�ӵ��ڽڵ�����������֮�����٣����ӷ���ʵ�����
%}

function [Sxy,AM,EdgeDelay,VertexDelay,ClusterMatrix]=cluster_topology(BorderLength,NodeAmount, ...
Alpha,Beta,VertexDegreeDUB,VertexSpeedDUB,VertexDirecDUB)
%% ��������б�
%BorderLength������������������ı߳�����λ��km
%NodeAmount������������ڵ�ĸ���
%Alpha����������������������AlphaԽ�󣬶̱���Գ��ߵı���Խ��
%Beta����������������������BetaԽ�󣬱ߵ��ܶ�Խ��
%PlotIf���������Ƿ���������ͼ�����Ϊ1����ͼ�����򲻻�ͼ

%%  �������
%Sxy��������3*N�ľ��󣬸��зֱ����ڴ洢�ڵ����ţ������꣬������ľ���
%AM��������0 1�洢����AM(i,j)=1��ʾ������i��j������ߣ�N*N
%EdgeDelay����������·ʱ�Ӿ���N*N
%VertexDelay���������ڵ�ʱ��������1*N
%%�Ƽ�������������� 
%BorderLength=1000;NodeAmount=25;Alpha=100000000;Beta=200000000000;
%PlotIf=1;
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

%���ڴ���·�������Ĳ���
condidatesNum = (VertexDegreeDUB(1) + VertexDegreeDUB(2))/4;

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
xlabel('x (km)','FontName','Times New Roman','FontSize',12);
ylabel('y (km)','FontName','Times New Roman','FontSize',12);

 %���������ʼ������Ԫ���飩
AM = cell(RowCnt,ColCnt);
EdgeDelay = cell(RowCnt,ColCnt);
VertexDelay = cell(RowCnt,ColCnt);
VertexMaxDegree = cell(RowCnt,ColCnt);
LDT = cell(RowCnt,ColCnt);
VertexStability = cell(RowCnt,ColCnt);
VertexPriority = cell(RowCnt,ColCnt);
IsClusterHead = cell(RowCnt,ColCnt);
LinkContri = cell(RowCnt,ColCnt);
Path = cell(RowCnt,ColCnt);

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


%������·�ĳ�ʼ����ֵ�����Ӿ���EdgeDelay)
for i = 1:RowCnt
    for j = 1:ColCnt
        Cluster = ClusterMatrix{i,j};
        nodesNum = size(Cluster,2);
        vertexMaxDegree = VertexMaxDegree{i,j};
        am = zeros(nodesNum,nodesNum);
        edgeDelay = Inf(nodesNum,nodesNum);
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
                    ldt(m,n) = GetLDT(MaxLinkDistance,Cluster(:,[m,n]));
                    ldt(n,m) = ldt(m,n);
                end
            end
            vertexStability(m) = sum(ldt(m,:));
        end
        vertexStability(nodesNum) = sum(ldt(nodesNum,:));
        
        AM{i,j} = am;
        EdgeDelay{i,j} = edgeDelay;
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
            vertexPriority(m) = GetPriority(m,vertexStability,am,vertexDelay);     
        end
        VertexPriority{i,j} = vertexPriority;
        IsClusterHead = SetClusterHead(IsClusterHead,[i,j],vertexPriority);
        Net_plot(BorderLength,nodesNum,ClusterMatrix{i,j},AM{i,j},...
            IsClusterHead{i,j},vertexPriority,[i,j]);
    end
end

%��ȡ�ڵ����·���׶�
for i = 1:RowCnt
    for j = 1:ColCnt
        am = AM{i,j};
        vertexMaxDegree = VertexMaxDegree{i,j};
        isClusterHead = IsClusterHead{i,j};
        edgeDelay = EdgeDelay{i,j};
        nodesNum = size(am,2);
        linkContri = zeros(1,nodesNum);
        for k = 1:nodesNum
            if isClusterHead(k) == 1
                headIdx = k;
                break;
            end
        end
        
        for k = 1:nodesNum
            if k ~= headIdx
                linkContri(k) = GetLinkContri(k,vertexMaxDegree(k),headIdx,am,...
                    edgeDelay);
            end
        end
        LinkContri{i,j} = linkContri;    
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
                [IsClusterHead,AM,EdgeDelay,VertexDelay,LDT,VertexStability,...
                    VertexPriority] = RemoveIneffectiveLink(IsClusterHead,...
                    MaxLinkDistance,k,ClusterMatrix,[i,j],AM,...
                    EdgeDelay,VertexDelay,LDT,VertexStability,VertexPriority);
                %���нڵ����룬��������
                [IsClusterHead,AM,EdgeDelay,VertexDelay,LDT,...
                    VertexStability,VertexPriority] ...
                    = ReConnect(IsClusterHead,MaxLinkDistance,ClusterMatrix,...
                    [i,j],k,AM,EdgeDelay,VertexDelay,...
                    VertexMaxDegree,LDT,VertexStability,VertexPriority);
                
                %�ڵ��ƶ���������(�ڵ�������˳���
                row = min(max(ceil(10 * Cluster(3,k)/(BorderLength*RowBase)),1),RowCnt);
                col = min(max(ceil(10 * Cluster(2,k)/(BorderLength*ColBase)),1),ColCnt);
                if (row ~= i) || (col ~= j)
                    [~,MostConnected] = CheckConnected(AM{i,j},k,inf);
                    %�ƶ��������صĽڵ㣬��������ͨ״̬���������
                    if MostConnected == 0
                        [IsClusterHead,ClusterMatrix,AM,EdgeDelay,...
                            VertexDelay,VertexMaxDegree,LDT,VertexStability,VertexPriority] ...
                            = VertexOutAndIn(IsClusterHead,MaxLinkDistance,ClusterMatrix,...
                            [i,j],k,Cluster(:,k),AM,EdgeDelay,VertexDelay,...
                            VertexMaxDegree,LDT,VertexStability,VertexPriority,...
                            BorderLength,RowBase,ColBase,RowCnt,ColCnt,VertexDegreeDUB);
                    end
                end
                
                Cluster = ClusterMatrix{i,j};
                nodesNum = size(Cluster,2);
                %���ص����һ��������룬��ʱk���
                if k > nodesNum
                    break;
                end
                
                %���벿�ֳ����븽��������
                [k,IsClusterHead,ClusterMatrix,AM,EdgeDelay,VertexDelay,VertexMaxDegree,...
                    LDT,VertexStability,VertexPriority]...
                    = ReConnectToOthers(IsClusterHead,MaxLinkDistance,ClusterMatrix,...
                    [i,j],k,AM,EdgeDelay,VertexDelay,VertexMaxDegree,...
                    LDT,VertexStability,VertexPriority,BorderLength,RowBase,ColBase,...
                    RowCnt,ColCnt);
                Cluster = ClusterMatrix{i,j};
                nodesNum = size(Cluster,2);
            end
            paths = ConstructPath(IsClusterHead{i,j},AM{i,j},LinkContri{i,j},...
                VertexStability{i,j},condidatesNum);
            
%             fprintf('======================Print Path======================\n')
%             fprintf('cluster[%d,%d]\n',i,j);
%             pathNum = size(paths,2);
%             for p = 1:pathNum
%                 path = paths{p};
%                 nodeNum = size(path,2);
%                 for q = 1:nodeNum
%                     fprintf('%d -> ',path(q));
%                 end
%                 fprintf('\n');
%             end
%             fprintf('\n==========================================\n');
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
            Net_plot(BorderLength,nodesNum,Cluster,AM{i,j},IsClusterHead{i,j},...
                VertexPriority{i,j},[i,j]);
        end
    end
end
    
end

%���ڻ����������˵ĺ���
function Net_plot(BorderLength,nodesNum,Cluster,am,isClusterHead,vertexPriority,...
                    ClusterIdx)
%���ڵ�
    plot(Cluster(2,:),Cluster(3,:),'ko','MarkerSize',5);    %'ko'����ɫԲȦ��'MarkerEdgeColor'����ǵı߿���ɫ��'MarkerFaceColor'����ǵ���ɫ
    hold on;    %��ʾ����ԭͼ���޸�
    %�ڵ�����
    for i = 1:nodesNum
        priorityStr = sprintf('%.3f',vertexPriority(i));
        Str = [int2str(i),'(',priorityStr,')'];
        if isClusterHead(i) == 1
            plot(Cluster(2,i),Cluster(3,i),'r.','MarkerSize',30);
        end
        text(Cluster(2,i)+BorderLength/100,Cluster(3,i)+BorderLength/100,Str,'FontName','Times New Roman','FontSize',12);
        hold on;
    end
%����
    lineColors = ['g','b','c','m','k'];
    num = size(lineColors,2);
    Idx = mod(2*ClusterIdx(1) + ClusterIdx(2),num); 
    if Idx == 0
        Idx = num;
    end
    lineColor = lineColors(Idx);
    
    for i = 1:(nodesNum-1)
        for j = (i+1):nodesNum
            if am(i,j) == 1
                plot([Cluster(2,i),Cluster(2,j)],[Cluster(3,i),Cluster(3,j)],lineColor);
                hold on;
            end
        end
    end
end

%�Ƴ���Ч��·
function [IsClusterHead,AM,EdgeDelay,VertexDelay,LDT,VertexStability,...
    VertexPriority] = RemoveIneffectiveLink(IsClusterHead_,MaxLinkDistance,nodeIdx,...
    ClusterMatrix_,ClusterIdx,AM_,EdgeDelay_,VertexDelay_,...
    LDT_,VertexStability_,VertexPriority_)

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
                [IsClusterHead_,AM_,EdgeDelay_,VertexDelay_,LDT_,...
                    VertexStability_,VertexPriority_] = LinkDisconnect(IsClusterHead_,...
                    [nodeIdx,i],ClusterIdx,AM_,EdgeDelay_,...
                    VertexDelay_,LDT_,VertexStability_,VertexPriority_);
            end
        end
    end
    IsClusterHead = IsClusterHead_;
    AM = AM_;
    EdgeDelay = EdgeDelay_;
    VertexDelay = VertexDelay_;
    LDT = LDT_;
    VertexStability = VertexStability_;
    VertexPriority = VertexPriority_;
end

%��·�Ͽ�
function [IsClusterHead,AM,EdgeDelay,VertexDelay,LDT,VertexStability,...
    VertexPriority] = LinkDisconnect(IsClusterHead_,LinkIdx,ClusterIdx,AM_,...
    EdgeDelay_,VertexDelay_,LDT_,VertexStability_,...
    VertexPriority_)

    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    edgeDelay = EdgeDelay_{ClusterIdx(1),ClusterIdx(2)};
    vertexDelay = VertexDelay_{ClusterIdx(1),ClusterIdx(2)};
    ldt = LDT_{ClusterIdx(1),ClusterIdx(2)};
    vertexStability = VertexStability_{ClusterIdx(1),ClusterIdx(2)};
    vertexPriority = VertexPriority_{ClusterIdx(1),ClusterIdx(2)};
    
    %����AM   
    am(LinkIdx(1),LinkIdx(2)) = 0;
    am(LinkIdx(2),LinkIdx(1)) = 0;
    AM_{ClusterIdx(1),ClusterIdx(2)} = am; 
    
    %����EdgeDelay
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
        vertexPriority(i) = GetPriority(i,vertexStability,am,vertexDelay);
    end
    VertexPriority_{ClusterIdx(1),ClusterIdx(2)} = vertexPriority;
    
    IsClusterHead = SetClusterHead(IsClusterHead_,ClusterIdx,vertexPriority);
    
    AM = AM_;
    EdgeDelay = EdgeDelay_;
    VertexDelay = VertexDelay_;
    LDT = LDT_;
    VertexStability = VertexStability_;
    VertexPriority = VertexPriority_;
end