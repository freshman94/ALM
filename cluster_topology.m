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


%���������������������ѡȡNN���ڵ�
for i = 1:NN
    SSxy(i,1) = BorderLength*rand;
    SSxy(i,2) = BorderLength*rand;
end
%�ڵ���ȷֲ�
[~,C] = kmeans(SSxy,NodeAmount);

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
Paths = cell(RowCnt,ColCnt);
InterClusterInfo = cell(RowCnt,ColCnt);
Tree = {};
Break_branches = {};
Src = 0;

%��ʼʱ��
t_pos = clock;
t_vd = clock;

%�ִ�(Cluster��������Ϣ�У���ţ��ڵ����꣬�ڵ��ٶȣ��ڵ㷽��,�ڵ������Ĵر��)
for i = 1:NodeAmount
    row = min(max(ceil(10 * Sxy(3,i)/(BorderLength*RowBase)),1),RowCnt);
    col = min(max(ceil(10 * Sxy(2,i)/(BorderLength * ColBase)),1),ColCnt);
    ClusterMatrix{row,col} = [ClusterMatrix{row,col},[Sxy(2,i),Sxy(3,i),Sxy(4,i),Sxy(5,i)]'];   
end

%Ϊ���ڽڵ��ţ�����ǽڵ������ر��
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
        clusterIdx = i+j;
        idxRow = ones(1,cols)*clusterIdx;
        tmp = [tmp;idxRow];
        ClusterMatrix{i,j} = tmp;
    end
end

%�ؼ以���ڵ����Ϣ���Ƿ����ӵ������أ����ӵ��Ĵر�ţ����ӵ��Ľڵ���ţ�
for i = 1:RowCnt
    for j = 1:ColCnt
        Cluster = ClusterMatrix{i,j};
        nodesNum = size(Cluster,2);
        infoNum = 3;
        info = cell(infoNum,nodesNum);
        for p = 1:nodesNum
           info{1,p} = 0; 
        end
        InterClusterInfo{i,j} = info;
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
        info = InterClusterInfo{i,j};
        nodesNum = size(vertexStability,2);
        vertexPriority = zeros(1,nodesNum);
        for m = 1:nodesNum
            vertexPriority(m) = GetPriority(m,vertexStability,am,vertexDelay,info);     
        end
        VertexPriority{i,j} = vertexPriority;
        IsClusterHead = SetClusterHead(IsClusterHead,[i,j],vertexPriority);
        Net_plot(BorderLength,nodesNum,ClusterMatrix,AM{i,j},...
            IsClusterHead{i,j},vertexPriority,[i,j],{},InterClusterInfo{i,j});
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
                    VertexPriority,InterClusterInfo] = RemoveIneffectiveLink...
                    (IsClusterHead,MaxLinkDistance,k,ClusterMatrix,[i,j],AM,EdgeDelay,...
                    VertexDelay,LDT,VertexStability,VertexPriority,InterClusterInfo);
                %���нڵ����룬��������
                [IsClusterHead,AM,EdgeDelay,VertexDelay,LDT,...
                    VertexStability,VertexPriority] ...
                    = ReConnect(IsClusterHead,MaxLinkDistance,ClusterMatrix,...
                    [i,j],k,AM,EdgeDelay,VertexDelay,VertexMaxDegree,LDT,...
                    VertexStability,VertexPriority,InterClusterInfo,RowCnt,ColCnt);
                
                %�ڵ��ƶ���������(�ڵ�������˳���
                row = min(max(ceil(10 * Cluster(3,k)/(BorderLength*RowBase)),1),RowCnt);
                col = min(max(ceil(10 * Cluster(2,k)/(BorderLength*ColBase)),1),ColCnt);
                if (row ~= i) || (col ~= j)
                    [~,MostConnected] = CheckConnected(AM{i,j},k,inf);
                    %�ƶ��������صĽڵ㣬��������ͨ״̬���������
                    if MostConnected == 0
                        [InterClusterInfo] = VertexOutAndIn(ClusterMatrix,...
                            MaxLinkDistance,[i,j],k,Cluster(:,k),AM,...
                            VertexMaxDegree,BorderLength,RowBase,ColBase,...
                            RowCnt,ColCnt,InterClusterInfo);
                    end
                end
                
                %�����븽��������
                [IsClusterHead,VertexPriority,InterClusterInfo]...
                    = ReConnectToOthers(IsClusterHead,MaxLinkDistance,ClusterMatrix,...
                    [i,j],k,AM,VertexDelay,VertexMaxDegree,...
                    VertexStability,VertexPriority,BorderLength,RowBase,...
                    ColBase,RowCnt,ColCnt,InterClusterInfo);
                k = k+1;
            end
            
            path = ConstructPath(IsClusterHead{i,j},AM{i,j},EdgeDelay{i,j});
            Paths{i,j} = path;
            
%             fprintf('======================Print Path======================\n')
%             fprintf('cluster[%d,%d]\n',i,j);
%             pathNum = size(path,2);
%             for p = 1:pathNum
%                 fprintf('%d ->%d\n ',path(1,p),path(2,p));
%             end
%             fprintf('\n==========================================\n');
        end
    end
    [Tree,Src,Break_branches] = ConstructMTree(ClusterMatrix,AM,InterClusterInfo,IsClusterHead,Paths,EdgeDelay);
    Src
    Tree
    cnt = size(Tree,1);
    for j = 1:cnt
       Tree{j} 
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
                
                %���ڵ㵽��صı߽磬�ڵ�Ӧ�����ķ����ƶ�
                %���ڵ����룬�ڵ�Ӧ�����ķ����ƶ�
                %(ʵ���У��ڵ��⵽����·�Ͽ���ʹ������ʱ����Ӧ�����ķ����ƶ�)
                BorderLeft = (j-1)*ColBase*BorderLength/10;
                BorderRight = j*ColBase*BorderLength/10;
                BorderTop = i*RowBase*BorderLength/10;
                BorderBtm = (i-1)*RowBase*BorderLength/10;
                [centerX,centerY] = GetCentroid(Cluster);
                am = AM{i,j};
                isClose = IsClose(Cluster);
                for k = 1:nodesNum
                    direc = atan2(centerY-Cluster(3,k),centerX-Cluster(2,k));
                    [~,MostConnected] = CheckConnected(am,k,inf);
                    if isClose == 1
                        Cluster(5,k) = direc + pi;
                    end
                    if Cluster(2,k) <= BorderLeft || Cluster(2,k) >= BorderRight...
                            || Cluster(3,k) <= BorderBtm || Cluster(3,k) >= BorderTop...
                            || MostConnected == 0                                      
                        Cluster(5,k) = direc;
                    end
                end
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
            Net_plot(BorderLength,nodesNum,ClusterMatrix,AM{i,j},IsClusterHead{i,j},...
                VertexPriority{i,j},[i,j],Paths{i,j},InterClusterInfo{i,j});
        end
    end
    Tree_plot(BorderLength,ClusterMatrix,Tree,Src,IsClusterHead);
    Break_plot(ClusterMatrix,Break_branches);
end
    
end

%���ڻ����������˵ĺ���
function Net_plot(BorderLength,nodesNum,ClusterMatrix,am,isClusterHead,vertexPriority,...
                    ClusterIdx,path,interCluserInfo)
    Colors = ['g','b','c','m'];
    num = size(Colors,2);
    Idx = mod(2*ClusterIdx(1) + ClusterIdx(2),num); 
    if Idx == 0
        Idx = num;
    end
    color = Colors(Idx);
    
    Cluster = ClusterMatrix{ClusterIdx(1),ClusterIdx(2)};
%���ڵ�
    nodeColor = strcat(color,'o');
    plot(Cluster(2,:),Cluster(3,:),nodeColor,'MarkerSize',5);    %'ko'����ɫԲȦ��'MarkerEdgeColor'����ǵı߿���ɫ��'MarkerFaceColor'����ǵ���ɫ
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
    for i = 1:(nodesNum-1)
        for j = (i+1):nodesNum
            if am(i,j) == 1
                plot([Cluster(2,i),Cluster(2,j)],[Cluster(3,i),Cluster(3,j)],'k');
                hold on;
            end
        end
    end
    
%���ؼ��
    for i = 1:nodesNum
        if interCluserInfo{1,i} == 1
            linkClusters = interCluserInfo{2,i};
            linkIdxs = interCluserInfo{3,i};
            linkNum = size(linkClusters,2);
            for j = 1:linkNum  ;
                linkClusterIdx = linkClusters(:,j);
                linkCluster = ClusterMatrix{linkClusterIdx(1),linkClusterIdx(2)};
                linkIdx = linkIdxs(j);
                plot([Cluster(2,i),linkCluster(2,linkIdx)],[Cluster(3,i),...
                    linkCluster(3,linkIdx)],'y');
            end
        end
    end
    %���ڴ���·��
    pathNum = size(path,2);
    for p = 1:pathNum
        if path(1,p) > 0
            plot([Cluster(2,path(1,p)),Cluster(2,path(2,p))],...
                [Cluster(3,path(1,p)),Cluster(3,path(2,p))],color);
%             [x,y] = GetArrow([Cluster(2,path(1,p)),Cluster(3,path(1,p))],...
%                 [Cluster(2,path(2,p)),Cluster(3,path(2,p))],BorderLength,[0.6,0.6]);
%             plot(x,y,color);
        end
    end
    
    hold on;     
end

%�����鲥��
function Tree_plot(BorderLength,ClusterMatrix,Tree,Src,IsClusterHead)
    isClusterHead = IsClusterHead{Src(1),Src(2)};
    nodesNum = size(isClusterHead,2);
    for i = 1:nodesNum
       if isClusterHead(i) == 1
           headIdx = i;
       end
    end
    cluster = ClusterMatrix{Src(1),Src(2)};
    Str = 'Src';
    text(cluster(2,headIdx),cluster(3,headIdx)+BorderLength/30,...
        Str,'FontName','Times New Roman','FontSize',12,'Color','red');
    
    num = size(Tree,1);
    for i = 1:2:(num-1)
       clusterIdx =  Tree{i};
       clusterNum = size(clusterIdx,1);
       path = Tree{i+1};
       %����
       if clusterNum == 1
           cluster = ClusterMatrix{clusterIdx(1),clusterIdx(2)};
           nodeNum = size(path,2);
           if nodeNum == 1
               continue;
           end
           for j = 1:(nodeNum-1)
               src = [cluster(2,path(j)),cluster(3,path(j))];
               dst = [cluster(2,path(j+1)),cluster(3,path(j+1))];
%                plot([src(1),dst(1)],[src(2),dst(2)],'r--');
               [x,y] = GetArrow(src,dst,BorderLength,[0.3,0.3]);
               plot(x,y,'r');
           end
       %�ؼ�
       elseif clusterNum == 2
           cluster_1 = ClusterMatrix{clusterIdx(1,1),clusterIdx(1,2)};
           cluster_2 = ClusterMatrix{clusterIdx(2,1),clusterIdx(2,2)};
           src = [cluster_1(2,path(1)),cluster_1(3,path(1))];
           dst = [cluster_2(2,path(2)),cluster_2(3,path(2))];
%            plot([src(1),dst(1)],[src(2),dst(2)],'r--');
           [x,y] = GetArrow(src,dst,BorderLength,[0.3,0.3]);
           plot(x,y,'r');
       end
    end
end

%�������벿�ִ���·��
function Break_plot(ClusterMatrix,Break_branches)    
    num = size(Break_branches,1);
    for i = 1:2:(num-1)
       clusterIdx =  Break_branches{i};
       clusterNum = size(clusterIdx,1);
       path = Break_branches{i+1};
       %����
       if clusterNum == 1
           cluster = ClusterMatrix{clusterIdx(1),clusterIdx(2)};
           nodeNum = size(path,2);
           if nodeNum == 1
               continue;
           end
           for j = 1:(nodeNum-1)
               src = [cluster(2,path(j)),cluster(3,path(j))];
               dst = [cluster(2,path(j+1)),cluster(3,path(j+1))];
               plot([src(1),dst(1)],[src(2),dst(2)],'y:');
%                [x,y] = GetArrow(src,dst,BorderLength,[0.3,0.3]);
%                plot(x,y,'r');
           end
       %�ؼ�
       elseif clusterNum == 2
           cluster_1 = ClusterMatrix{clusterIdx(1,1),clusterIdx(1,2)};
           cluster_2 = ClusterMatrix{clusterIdx(2,1),clusterIdx(2,2)};
           src = [cluster_1(2,path(1)),cluster_1(3,path(1))];
           dst = [cluster_2(2,path(2)),cluster_2(3,path(2))];
           plot([src(1),dst(1)],[src(2),dst(2)],'y:');
%            [x,y] = GetArrow(src,dst,BorderLength,[0.3,0.3]);
%            plot(x,y,'r');
       end
    end
end

%�Ƴ���Ч��·
function [IsClusterHead,AM,EdgeDelay,VertexDelay,LDT,VertexStability,VertexPriority,...
    InterClusterInfo] = RemoveIneffectiveLink(IsClusterHead_,MaxLinkDistance,nodeIdx,...
    ClusterMatrix_,ClusterIdx,AM_,EdgeDelay_,VertexDelay_,...
    LDT_,VertexStability_,VertexPriority_,InterClusterInfo_)

    Cluster = ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)};
    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    
    %���Բ���
    once = 1;
    %�Ƴ�������Ч��·
    nodesNum = size(am,2);
    for i = 1:nodesNum
        if am(nodeIdx,i) == 1 
            distance = ( (Cluster(2,nodeIdx)-Cluster(2,i))^2 +...
                    (Cluster(3,nodeIdx)-Cluster(3,i))^2)^0.5;
            if distance > MaxLinkDistance
                %����
                if once == 1
%                     fprintf('================================RemoveIneffectiveLink==================\n');
%                     fprintf('clusterIdx=[%d,%d]\n',ClusterIdx(1),ClusterIdx(2));
%                     for m = 1:nodesNum
%                         for n = 1:nodesNum
%                             fprintf('am(%d,%d) = %d\t',m,n,am(m,n));
%                         end
%                         fprintf('\n');
%                     end
                end
                once = 2;
%                 fprintf('RemoveLink (%d,%d)\n',nodeIdx,i);
                [IsClusterHead_,AM_,EdgeDelay_,VertexDelay_,LDT_,...
                    VertexStability_,VertexPriority_] = InnerLinkDisconnect(IsClusterHead_,...
                    [nodeIdx,i],ClusterIdx,AM_,EdgeDelay_,...
                    VertexDelay_,LDT_,VertexStability_,VertexPriority_,InterClusterInfo_);
            end
        end
    end
    
    %�Ƴ��ؼ���Ч��·
    [InterClusterInfo_,VertexPriority_,IsClusterHead_] = InterLinkDisconnect...
            (InterClusterInfo_,VertexPriority_,ClusterMatrix_,ClusterIdx,nodeIdx,...
            MaxLinkDistance,IsClusterHead_,AM_,VertexDelay_,VertexStability_);
    
    
    IsClusterHead = IsClusterHead_;
    AM = AM_;
    EdgeDelay = EdgeDelay_;
    VertexDelay = VertexDelay_;
    LDT = LDT_;
    VertexStability = VertexStability_;
    VertexPriority = VertexPriority_;
    InterClusterInfo = InterClusterInfo_;
end

%��·�Ͽ�
function [IsClusterHead,AM,EdgeDelay,VertexDelay,LDT,VertexStability,...
    VertexPriority] = InnerLinkDisconnect(IsClusterHead_,LinkIdx,ClusterIdx,AM_,...
    EdgeDelay_,VertexDelay_,LDT_,VertexStability_,...
    VertexPriority_,InterClusterInfo_)

    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    edgeDelay = EdgeDelay_{ClusterIdx(1),ClusterIdx(2)};
    vertexDelay = VertexDelay_{ClusterIdx(1),ClusterIdx(2)};
    ldt = LDT_{ClusterIdx(1),ClusterIdx(2)};
    vertexStability = VertexStability_{ClusterIdx(1),ClusterIdx(2)};
    vertexPriority = VertexPriority_{ClusterIdx(1),ClusterIdx(2)};
    info = InterClusterInfo_{ClusterIdx(1),ClusterIdx(2)};
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
        vertexPriority(i) = GetPriority(i,vertexStability,am,vertexDelay,info);
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

function [InterClusterInfo,VertexPriority,IsClusterHead] = InterLinkDisconnect...
            (InterClusterInfo_,VertexPriority_,ClusterMatrix,ClusterIdx,nodeIdx,...
            MaxLinkDistance,IsClusterHead_,AM_,VertexDelay_,VertexStability_)
        info = InterClusterInfo_{ClusterIdx(1),ClusterIdx(2)};
        if info{1,nodeIdx} == 1
            linkClusters = info{2,nodeIdx};
            linkNum = size(linkClusters,2);
            linkIdxs = info{3,nodeIdx};
            Cluster = ClusterMatrix{ClusterIdx(1),ClusterIdx(2)};
            nodePos = Cluster(:,nodeIdx);
            vertexPriority = VertexPriority_{ClusterIdx(1),ClusterIdx(2)};
            am = AM_{ClusterIdx(1),ClusterIdx(2)};
            vertexDelay = VertexDelay_{ClusterIdx(1),ClusterIdx(2)};
            vertexStability = VertexStability_{ClusterIdx(1),ClusterIdx(2)};
            nodesNum = size(am,2);
            i = 1;
            while i <= linkNum
                linkClusterIdx = linkClusters(:,i);
                linkIdx = linkIdxs(i);
                linkCluster = ClusterMatrix{linkClusterIdx(1),linkClusterIdx(2)};
                linkNodePos = linkCluster(:,linkIdx);
                Distance = ( (nodePos(2)-linkNodePos(2))^2 +...
                    (nodePos(3)-linkNodePos(3))^2)^0.5;
                if Distance > MaxLinkDistance
%                     fprintf('Remove inter link([%d,%d]%d,[%d,%d]%d) \n', ClusterIdx(1),...
%                         ClusterIdx(2),nodeIdx,linkClusterIdx(1),linkClusterIdx(2),linkIdx);
%                     info
%                     nodeIdx
%                     linkClusters
%                     linkIdxs
                    %���½ڵ�nodeIdx��Ϣ
                    linkClusters(:,i) = [];
                    linkIdxs(i) = [];
                    info{2,nodeIdx} = linkClusters;
                    info{3,nodeIdx} = linkIdxs;
                    for k = 1:nodesNum
                        vertexPriority(k) = GetPriority(k,vertexStability,am,...
                            vertexDelay,info);
                    end
                    IsClusterHead_ = SetClusterHead(IsClusterHead_,ClusterIdx,vertexPriority);
                    
                    %���½��linkIdx��Ϣ
                    vertexPriority_ = VertexPriority_{linkClusterIdx(1),linkClusterIdx(2)};
                    am_ = AM_{linkClusterIdx(1),linkClusterIdx(2)};
                    vertexDelay_ = VertexDelay_{linkClusterIdx(1),linkClusterIdx(2)};
                    vertexStability_ = VertexStability_{linkClusterIdx(1),linkClusterIdx(2)};
                    linkInfo = InterClusterInfo_{linkClusterIdx(1),linkClusterIdx(2)};
                    linkInfo_clusters = linkInfo{2,linkIdx};
                    linkInfo_idxs = linkInfo{3,linkIdx};
                    num = size(linkInfo_idxs,2);
                    k = 1;
                    while k <= num
                        if linkInfo_clusters(1,k) == ClusterIdx(1) && linkInfo_clusters(2,k)...
                            == ClusterIdx(2) && linkInfo_idxs(k) == nodeIdx
%                             linkInfo
%                             linkInfo_clusters
%                             linkInfo_idxs
                            linkInfo_clusters(:,k) = [];
                            linkInfo_idxs(k) = [];
                            linkInfo{2,linkIdx} = linkInfo_clusters;
                            linkInfo{3,linkIdx} = linkInfo_idxs;
                            num = size(linkInfo_idxs,2);
                        else
                            k = k+1;
                        end
                    end
                    for k = 1:num
                        vertexPriority_(k) = GetPriority(k,vertexStability_,am_,...
                            vertexDelay_,linkInfo);
                    end
                    IsClusterHead_ = SetClusterHead(IsClusterHead_,linkClusterIdx,vertexPriority_);
                    if size(linkInfo{2,linkIdx},2) == 0
                        linkInfo{1,linkIdx} = 0;
                    end
                    InterClusterInfo_{linkClusterIdx(1),linkClusterIdx(2)} = linkInfo;
                    VertexPriority_{linkClusterIdx(1),linkClusterIdx(2)} = vertexPriority_;
                
                    linkNum = size(linkClusters,2);
                    if size(info{2,nodeIdx},2) == 0
                        info{1,nodeIdx} = 0;
                        break;
                    end
                else
                    i = i+1;
                end
                
            end
            InterClusterInfo_{ClusterIdx(1),ClusterIdx(2)} = info;
            VertexPriority_{ClusterIdx(1),ClusterIdx(2)} = vertexPriority;
        end
        InterClusterInfo = InterClusterInfo_;
        VertexPriority = VertexPriority_;    
        IsClusterHead = IsClusterHead_;
end

function [centerX,centerY] = GetCentroid(Cluster)
    nodesNum = size(Cluster,2);
    Xsum = 0;
    Ysum = 0;
    for i = 1:nodesNum
        Xsum = Xsum + Cluster(2,i);
        Ysum = Ysum + Cluster(3,i);
    end
    centerX = Xsum/nodesNum;
    centerY = Ysum/nodesNum;
end

function [isClose] = IsClose(Cluster)
    [centerX,centerY] = GetCentroid(Cluster);
    nodesNum = size(Cluster,2);
    variance = 0;
    thr = 2200;
    isClose = 0;
    for i = 1:nodesNum
        variance = variance+(Cluster(2,i)-centerX).^2 + (Cluster(3,i)-centerY).^2;
    end
    if variance < thr
        isClose = 1;
    end
end