function [IsClusterHead,ClusterMatrix,AM,EdgeDelay,VertexDelay,...
    VertexMaxDegree,LDT,VertexStability,VertexPriority]...
    = VertexOutAndIn(IsClusterHead_,MaxLinkDistance,ClusterMatrix_,ClusterIdx,nodeIdx,...
    nodePos,AM_,EdgeDelay_,VertexDelay_,VertexMaxDegree_,LDT_,VertexStability_,...
    VertexPriority_,BorderLength,RowBase,ColBase,RowCnt,ColCnt,VertexDegreeDUB)
%% ���������Ϊ��
%% ��1��������Ӱ��ڵ����Ϣ
%% ��2����Ӱ��ڵ㳢������
%% ��3���ڵ��뿪��ɾ���ڵ����Ϣ��
%% ��4���ڵ������һ����

%%�������
%ClusterIdx: 1*2��������ָ��ĳ����
%nodeIdx: ��ֵ��ָ������ĳ���ڵ�

    ClusterOut = ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)};
    nodesOutNum = size(ClusterOut,2);
    %ͳ����Ӱ��ڵ�
    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    %�ж���Ӱ��ڵ��Ƿ���Ҫ���¼����
    affectedNodes = [];
    for i = 1:nodesOutNum
        if am(nodeIdx,i) == 1 
            affectedNodes = [affectedNodes,i];
        end
    end  
    
    affectedNum = size(affectedNodes,2);
    
    %����
    fprintf('clusterIdx=[%d,%d],nodeIdx = %d,affectedNum = %d\n',ClusterIdx(1),ClusterIdx(2),nodeIdx,affectedNum);
    
    fprintf('================print am==========\n');
    for i = 1:nodesOutNum
        for j = 1:nodesOutNum
            fprintf('am(%d,%d) = %d\t',i,j,am(i,j));
        end
        fprintf('\n');
    end
    fprintf('=========================\n');
    
    %������Ӱ������Ϣ
    [AM_,LDT_,VertexDelay_,VertexStability_,VertexPriority_,IsClusterHead_] = ...
    UpdateInfo(ClusterIdx,nodeIdx,affectedNodes,AM_,LDT_,VertexDelay_,EdgeDelay_,...
    VertexStability_,VertexMaxDegree_,VertexPriority_,IsClusterHead_);
    
    %��Ӱ��ڵ㳢������
    %�˱�������ָʾ����Ľڵ��Ƿ����¼���ɹ�
    AllConnected = 1;
    for i = 1:affectedNum
        am = AM_{ClusterIdx(1),ClusterIdx(2)};
        %���ڵ�Ķ�Ϊ1�������˳���������������Ӱ��
        isConnected = CheckConnected(am,affectedNodes(i),nodeIdx);
        %���Դ���
        fprintf('affectedNode = %d,isConnected=%d\n',affectedNodes(i),isConnected);
        %�ýڵ����룬���¼���
        if isConnected == 0
            affectedNodePos = ClusterOut(:,affectedNodes(i));
            [IsClusterHead_,AllConnected,AM_,EdgeDelay_,VertexDelay_,LDT_,...
                VertexStability_,VertexPriority_] = ...
            OldNodeJoin(IsClusterHead_,MaxLinkDistance,affectedNodes,ClusterMatrix_,...
            ClusterIdx,affectedNodes(i),affectedNodePos,nodeIdx,AM_,...
            EdgeDelay_,VertexDelay_,VertexMaxDegree_,LDT_,VertexStability_,...
            VertexPriority_);
            %������ͨ���ؽ��ɹ�
%             if AllConnected == 1
%                 break;
%             end
        end
    end
    
    fprintf('AllConnected = %d\n',AllConnected);
    
    %����Ľڵ�û�м���ɹ�����nodeIdx�������ڽڵ������������¹���������ͨ��
    if AllConnected == 0
        [IsClusterHead_,AM_,EdgeDelay_,VertexDelay_,LDT_,...
        VertexStability_,VertexPriority_] = ...
        ConnectNodeIdxAllAdj(IsClusterHead_,MaxLinkDistance,ClusterMatrix_,...
        ClusterIdx,affectedNodes,AM_,EdgeDelay_,VertexDelay_,...
        VertexMaxDegree_,LDT_,VertexStability_,VertexPriority_);
    end
    IsClusterHead = IsClusterHead_;
    AM = AM_;
    EdgeDelay = EdgeDelay_;
    VertexDelay = VertexDelay_;
    LDT = LDT_;
    VertexStability = VertexStability_;
    VertexPriority = VertexPriority_;
    
    %�ڵ��뿪(ɾ���ڵ���Ϣ)
    [IsClusterHead_,ClusterMatrix_,AM_,EdgeDelay_,VertexDelay_,...
        VertexMaxDegree_,LDT_,VertexStability_,VertexPriority_]...
    = NodeOut(IsClusterHead_,nodeIdx,nodesOutNum,ClusterMatrix_,ClusterIdx,AM_,...
        EdgeDelay_,VertexDelay_,VertexMaxDegree_,LDT_,...
        VertexStability_,VertexPriority_);
    
    %%�ڵ����
    [IsClusterHead,ClusterMatrix,AM,EdgeDelay,VertexDelay,...
        VertexMaxDegree,LDT,VertexStability,VertexPriority] = ...
    NewNodeJoin(IsClusterHead_,MaxLinkDistance,ClusterMatrix_,nodePos,AM_,...
    EdgeDelay_,VertexDelay_,VertexMaxDegree_,LDT_,VertexStability_,...
    VertexPriority_,BorderLength,RowBase,ColBase,RowCnt,ColCnt,VertexDegreeDUB);
end

%������Ӱ��ڵ����Ϣ
function [AM,LDT,VertexDelay,VertexStability,VertexPriority,IsClusterHead] = ...
    UpdateInfo(ClusterIdx,nodeIdx,affectedNodes,AM_,LDT_,VertexDelay_,EdgeDelay,...
    VertexStability_,VertexMaxdegree,VertexPriority_,IsClusterHead_)
    
    fprintf('===============================UpdateInfo======================\n');
    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    ldt = LDT_{ClusterIdx(1),ClusterIdx(2)};
    vertexDelay = VertexDelay_{ClusterIdx(1),ClusterIdx(2)};
    edgeDelay = EdgeDelay{ClusterIdx(1),ClusterIdx(2)};
    vertexStability = VertexStability_{ClusterIdx(1),ClusterIdx(2)};
    vertexMaxDegree = VertexMaxdegree{ClusterIdx(1),ClusterIdx(2)};
    vertexPriority = VertexPriority_{ClusterIdx(1),ClusterIdx(2)};
    nodesNum = size(am,1);
    affectedNodesNum = size(affectedNodes,2);
    
    am(:,nodeIdx) = zeros(nodesNum,1);
    am(nodeIdx,:) = zeros(1,nodesNum);
    
    ldt(:,nodeIdx) = zeros(nodesNum,1);
    ldt(nodeIdx,:) = zeros(1,nodesNum);
    vertexDelay(nodeIdx) = 0;
    for i = 1:affectedNodesNum
        node = affectedNodes(i);
        vertexStability(node) = sum(ldt(node,:));
    end
    
    for i = 1:nodesNum
        vertexDelay(i) = GetVertexDelay(i,am,edgeDelay);
    end
    
    for i = 1:nodesNum
        vertexPriority(i) = GetPriority(i,vertexStability,am,vertexDelay);
    end
    
    IsClusterHead = SetClusterHead(IsClusterHead_,ClusterIdx,vertexPriority);
    
    AM_{ClusterIdx(1),ClusterIdx(2)} = am;
    LDT_{ClusterIdx(1),ClusterIdx(2)} = ldt;
    VertexDelay_{ClusterIdx(1),ClusterIdx(2)} = vertexDelay;
    VertexStability_{ClusterIdx(1),ClusterIdx(2)} = vertexStability;
    VertexPriority_{ClusterIdx(1),ClusterIdx(2)} = vertexPriority;
    
    AM = AM_;
    LDT = LDT_;
    VertexDelay = VertexDelay_;
    VertexStability = VertexStability_;
    VertexPriority = VertexPriority_;  
end

%%��Ӱ��ڵ㳢������
%ͬһ���ؽڵ����¼���
function [IsClusterHead,Connected,AM,EdgeDelay,VertexDelay,LDT,...
    VertexStability,VertexPriority] = ...
    OldNodeJoin(IsClusterHead_,MaxLinkDistance,affectedNodes,ClusterMatrix_,ClusterIdx,...
    affectedNodeIdx,affectedNodePos,nodeIdx,AM_,EdgeDelay_,...
    VertexDelay_,VertexMaxDegree_,LDT_,VertexStability_,VertexPriority_)
%% ���˳��ڵ�����ڽڵ���ѡ���������Ľڵ�����

    fprintf('===========================OldNodeJoin===========================\n');
    %ѡ����nodeIdx������ڽ������
    sortedNodes = SortNodesByDistance(0,affectedNodePos,ClusterMatrix_,ClusterIdx,affectedNodes);
    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    vertexMaxDegree = VertexMaxDegree_{ClusterIdx(1),ClusterIdx(2)};
    affectedNum = size(sortedNodes,2);
    Connected = 0;
    Cluster = ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)};  
    for i = 1:affectedNum
        linkIdx = sortedNodes(1,i);
        linkNodeDegree = sum(am(linkIdx,:));
        linkNodeMaxDegree = vertexMaxDegree(linkIdx);
        if linkIdx ~= affectedNodeIdx && am(affectedNodeIdx,linkIdx) == 0
            if linkNodeDegree >= linkNodeMaxDegree
                continue;
            end
            Distance = sortedNodes(2,i);
            if Distance > MaxLinkDistance
                break;
            end
            
            am(affectedNodeIdx,linkIdx) = 1;
            am(linkIdx,affectedNodeIdx) = 1;
            
            %����
            fprintf('linkIdx = %d\n',linkIdx);
            nodes = size(am,2);
            for k = 1:nodes
                fprintf('am(%d,%d) = %d\n',affectedNodeIdx,k,am(affectedNodeIdx,k));
            end
            
            Connected = CheckConnected(am,affectedNodeIdx,nodeIdx);
            fprintf('Connected=%d\n',Connected);
            if Connected == 1
                break;
            else
                am(affectedNodeIdx,linkIdx) = 0;
                am(linkIdx,affectedNodeIdx) = 0;
            end
        end
    end
    
    %�ҵ����ӵĽ��ʹ������ͨ
    if Connected == 1     
        %����AM
        AM_{ClusterIdx(1),ClusterIdx(2)} = am;

        %����EdgeDelay
        edgeDelay = EdgeDelay_{ClusterIdx(1),ClusterIdx(2)};
        NewEdgeDelay = 0.5*Distance/100000;
        edgeDelay(affectedNodeIdx,linkIdx) = NewEdgeDelay;
        edgeDelay(linkIdx,affectedNodeIdx) = edgeDelay(affectedNodeIdx,linkIdx);
        EdgeDelay_{ClusterIdx(1),ClusterIdx(2)} = edgeDelay;

        %����VertexDelay
        vertexDelay = VertexDelay_{ClusterIdx(1),ClusterIdx(2)};
        nodesNum = size(am,2);
        for i = 1:nodesNum
            vertexDelay(i) = GetVertexDelay(i,am,edgeDelay);
        end
        VertexDelay_{ClusterIdx(1),ClusterIdx(2)} = vertexDelay;
        
        %����LDT
        ldt = LDT_{ClusterIdx(1),ClusterIdx(2)};
        ldt(affectedNodeIdx,linkIdx) = GetLDT(MaxLinkDistance,...
                                        Cluster(:,[affectedNodeIdx,linkIdx]));
        ldt(linkIdx,affectedNodeIdx) = ldt(affectedNodeIdx,linkIdx);
        LDT_{ClusterIdx(1),ClusterIdx(2)} = ldt;
        
        %����VertexStability
        vertexStability = VertexStability_{ClusterIdx(1),ClusterIdx(2)};
        vertexStability(affectedNodeIdx) = sum(ldt(affectedNodeIdx,:));
        vertexStability(linkIdx) = sum(ldt(linkIdx,:));
        VertexStability_{ClusterIdx(1),ClusterIdx(2)} = vertexStability;
        
        vertexPriority = VertexPriority_{ClusterIdx(1),ClusterIdx(2)};
        for i = 1:nodesNum
            vertexPriority(i) = GetPriority(i,vertexStability,am,vertexDelay);
        end
        VertexPriority_{ClusterIdx(1),ClusterIdx(2)} = vertexPriority;
        IsClusterHead_ = SetClusterHead(IsClusterHead_,ClusterIdx,vertexPriority);   
    end
    AM = AM_;
    EdgeDelay = EdgeDelay_;
    VertexDelay = VertexDelay_;
    LDT = LDT_;
    VertexStability = VertexStability_;
    VertexPriority = VertexPriority_;
    IsClusterHead = IsClusterHead_;
end

function [IsClusterHead,AM,EdgeDelay,VertexDelay,LDT,VertexStability,VertexPriority]...
    = ConnectNodeIdxAllAdj(IsClusterHead_,MaxLinkDistance,...
    ClusterMatrix_,ClusterIdx,affectedNodes,AM_,EdgeDelay_,...
    VertexDelay_,VertexMaxDegree_,LDT_,VertexStability_,VertexPriority_)
%%��affectedNodes�������������Ȼ����β����
    fprintf('=========================ConnectNodeIdxAllAdj===================================\n');
    affectedNum = size(affectedNodes,2);
    Cluster = ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)};
    
    affectedNodesWithX = zeros(2,affectedNum);
    for i = 1:affectedNum
        affectedNodesWithX(:,i) = Cluster(1:2,affectedNodes(i));
    end
    X = affectedNodesWithX(2,:);
    sorted_X = sort(X);
    sortedAffectedNodes = zeros(1,affectedNum);
    for i = 1:affectedNum
        pos = find(X==sorted_X(i));
        sortedAffectedNodes(i) = affectedNodesWithX(1,pos);
    end
    
    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    vertexMaxDegree = VertexMaxDegree_{ClusterIdx(1),ClusterIdx(2)};
    edgeDelay = EdgeDelay_{ClusterIdx(1),ClusterIdx(2)};
    ldt = LDT_{ClusterIdx(1),ClusterIdx(2)};
    vertexDelay = VertexDelay_{ClusterIdx(1),ClusterIdx(2)};
    vertexStability = VertexStability_{ClusterIdx(1),ClusterIdx(2)};
    vertexPriority = VertexPriority_{ClusterIdx(1),ClusterIdx(2)};
    
    node = sortedAffectedNodes(1);
    for i = 2:affectedNum
        %ֻ���������Ľ��(��·��������Ч,�����Ƕ�Լ����
        nodeDegree = sum(am(node,:));
        nodeMaxDegree = vertexMaxDegree(node);
        nodePos = Cluster(:,node);
        nextNode = sortedAffectedNodes(i);
        nextNodeDegree = sum(am(nextNode,:));
        nextNodeMaxDegree = vertexMaxDegree(nextNode);
        nextNodePos = Cluster(:,nextNode);
        if nodeDegree >= nodeMaxDegree || ...
            nextNodeDegree >= nextNodeMaxDegree
            node = nextNode;
            continue;
        end
        
        Distance = ( (nodePos(2)-nextNodePos(2))^2 +...
                        (nodePos(3)-nextNodePos(3))^2)^0.5;
        if am(node,nextNode) == 0 && Distance < MaxLinkDistance   
            %����AM_
            am(node,nextNode) = 1;
            am(nextNode,node) = 1;      

            %����EdgeDelay
            NewEdgeDelay = 0.5*Distance/100000;
            edgeDelay(node,nextNode) = NewEdgeDelay;
            edgeDelay(nextNode,node) = edgeDelay(node,nextNode);
                      
            %����VertexDelay
            nodesNum = size(am,2);
            for k = 1:nodesNum
                vertexDelay(k) = GetVertexDelay(k,am,edgeDelay);
            end         
            
            %����LDT
            ldt(node,nextNode) = GetLDT(MaxLinkDistance,Cluster(:,[node,nextNode]));
            ldt(nextNode,node) = ldt(node,nextNode);
            
            %����VertexStability
            vertexStability(node) = sum(ldt(node,:));
            vertexStability(nextNode) = sum(ldt(nextNode,:));
            
            %����VertexPriority
            for k = 1:nodesNum
                vertexPriority(i) = GetPriority(i,vertexStability,am,vertexDelay);
            end
            IsClusterHead_ = SetClusterHead(IsClusterHead_,ClusterIdx,vertexPriority);
            
        end
        node = nextNode;
    end
    
    %����
    fprintf('===================================PrintAM====================\n');
    nodes = size(am,2);
    for i = 1:nodes
        for j = 1:nodes
            fprintf('am(%d,%d) = %d\t',i,j,am(i,j));
        end
        fprintf('\n');
    end
    fprintf('==============================================\n');
    
    AM_{ClusterIdx(1),ClusterIdx(2)} = am;
    EdgeDelay_{ClusterIdx(1),ClusterIdx(2)} = edgeDelay;
    LDT_{ClusterIdx(1),ClusterIdx(2)} = ldt;
    VertexDelay_{ClusterIdx(1),ClusterIdx(2)} = vertexDelay; 
    VertexStability_{ClusterIdx(1),ClusterIdx(2)} = vertexStability;
    VertexPriority_{ClusterIdx(1),ClusterIdx(2)} = vertexPriority;
    
    AM = AM_;
    EdgeDelay = EdgeDelay_;
    VertexDelay = VertexDelay_;
    LDT = LDT_;
    VertexStability = VertexStability_;
    VertexPriority = VertexPriority_;
    IsClusterHead = IsClusterHead_;
end

%�ڵ��뿪
function [IsClusterHead,ClusterMatrix,AM,EdgeDelay,VertexDelay,...
    VertexMaxDegree,LDT,VertexStability,VertexPriority]...
    = NodeOut(IsClusterHead_,nodeOutIdx,nodesOutNum,ClusterMatrix_,ClusterIdx,AM_,...
    EdgeDelay_,VertexDelay_,VertexMaxDegree_,LDT_,VertexStability_,...
    VertexPriority_)
 
    fprintf('===========================NodeOut===========================\n');
    fprintf('nodeOutIdx = %d\n',nodeOutIdx);
    Cluster = ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)};
    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    edgeDelay = EdgeDelay_{ClusterIdx(1),ClusterIdx(2)};
    vertexDelay = VertexDelay_{ClusterIdx(1),ClusterIdx(2)};
    vertexMaxDegree = VertexMaxDegree_{ClusterIdx(1),ClusterIdx(2)};
    ldt = LDT_{ClusterIdx(1),ClusterIdx(2)};
    vertexStability = VertexStability_{ClusterIdx(1),ClusterIdx(2)};
    vertexPriority = VertexPriority_{ClusterIdx(1),ClusterIdx(2)};
    isClusterHead = IsClusterHead_{ClusterIdx(1),ClusterIdx(2)};
    
    %����AM   
    am(:,nodeOutIdx) = [];
    am(nodeOutIdx,:) = [];
    AM_{ClusterIdx(1),ClusterIdx(2)} = am; 
    
    %����ClusterMatrix
    Cluster(:,nodeOutIdx) = [];
    Cluster(1,:) = [1:nodesOutNum-1];
    ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)} = Cluster;
    
    
    %����EdgeDelay
    edgeDelay(:,nodeOutIdx) = [];
    edgeDelay(nodeOutIdx,:) = [];
    EdgeDelay_{ClusterIdx(1),ClusterIdx(2)} = edgeDelay; 
    
    %����LDT
    ldt(:,nodeOutIdx) = [];
    ldt(nodeOutIdx,:) = [];
    LDT_{ClusterIdx(1),ClusterIdx(2)} = ldt;
    
    %����Vertex*    
    vertexDelay(nodeOutIdx) = [];
    VertexDelay_{ClusterIdx(1),ClusterIdx(2)} = vertexDelay;
    
    vertexMaxDegree(nodeOutIdx) = [];
    VertexMaxDegree_{ClusterIdx(1),ClusterIdx(2)} = vertexMaxDegree;
    
    vertexStability(nodeOutIdx) = [];
    VertexStability_{ClusterIdx(1),ClusterIdx(2)} = vertexStability;
    
    vertexPriority(nodeOutIdx) = [];
    VertexPriority_{ClusterIdx(1),ClusterIdx(2)} = vertexPriority;
    
    %�뿪�Ľڵ�Ϊ����
    if isClusterHead(nodeOutIdx) == 1
        isClusterHead(nodeOutIdx) = [];
        IsClusterHead_{ClusterIdx(1),ClusterIdx(2)} = isClusterHead;
        IsClusterHead_ = SetClusterHead(IsClusterHead_,ClusterIdx,vertexPriority);
    else
        isClusterHead(nodeOutIdx) = [];
        IsClusterHead_{ClusterIdx(1),ClusterIdx(2)} = isClusterHead;
    end
    
    ClusterMatrix = ClusterMatrix_;
    AM = AM_;
    EdgeDelay = EdgeDelay_;
    VertexDelay = VertexDelay_;
    VertexMaxDegree = VertexMaxDegree_;
    LDT = LDT_;
    VertexStability = VertexStability_;
    VertexPriority = VertexPriority_;
    IsClusterHead = IsClusterHead_;
end


%��һ���ؽڵ����
function [IsClusterHead,ClusterMatrix,AM,EdgeDelay,VertexDelay,...
    VertexMaxDegree,LDT,VertexStability,VertexPriority]...
    = NewNodeJoin(IsClusterHead_,MaxLinkDistance,ClusterMatrix_,nodePos,AM_,...
    EdgeDelay_,VertexDelay_,VertexMaxDegree_,LDT_,VertexStability_,...
    VertexPriority_,BorderLength,RowBase,ColBase,RowCnt,ColCnt,VertexDegreeDUB)
    
    fprintf('===========================NewNodeJoin===========================\n');
    row = min(max(ceil(10 * nodePos(3)/(BorderLength*RowBase)),1),RowCnt);
    col = min(max(ceil(10 * nodePos(2)/(BorderLength * ColBase)),1),ColCnt);

    %����ClusterMatrix
    Cluster = ClusterMatrix_{row,col};
    nodesNum = size(Cluster,2);
    Cluster = [Cluster,[nodesNum+1,nodePos(2),nodePos(3),nodePos(4),nodePos(5)]'];
    ClusterMatrix_{row,col} = Cluster;
    ClusterMatrix = ClusterMatrix_;

    am = AM_{row,col};
    amAddCol = zeros(nodesNum,1);
    amAddRow = zeros(1,nodesNum+1);
    
    edgeDelay = EdgeDelay_{row,col};
    edgeDelayAddCol = Inf(nodesNum,1);
    edgeDelayAddRow = Inf(1,nodesNum+1);
    
    vertexDelay = VertexDelay_{row,col};
    vertexDelay(nodesNum+1) = 10e10;
    
    vertexMaxDegree = VertexMaxDegree_{row,col};
    vertexMaxDegree(nodesNum+1) = VertexDegreeDUB(1)+round( (VertexDegreeDUB(2)-VertexDegreeDUB(1))*rand);
    VertexMaxDegree_{row,col} = vertexMaxDegree;
    
    ldt = LDT_{row,col};
    ldtAddCol = zeros(nodesNum,1);
    ldtAddRow = zeros(1,nodesNum+1);
    
    vertexStability = VertexStability_{row,col};
    vertexStability(nodesNum+1) = 0;
    
    vertexPriority = VertexPriority_{row,col};
    vertexPriority(nodesNum+1) = 1;
    %ѡ���������Ľڵ���·(�����·��Ч�ԣ�
    %��1�������linkIdx�Ķȳ������ֵ����ѡ����һ��㳢������
    % (2��������ľ��볬����δ����������Ӧ�����ò���ֵ���������break
    sortedNodes = SortNodesByDistance(1,nodePos,ClusterMatrix_,[row,col],[1:nodesNum]);
    for i = 1:nodesNum  
        linkIdx = sortedNodes(1,i);
        linkNodePos = Cluster(:,linkIdx);
        linkNodeDegree = sum(am(linkIdx,:));
        linkNodeMaxDegree = vertexMaxDegree(linkIdx);
        Distance = sortedNodes(2,i);
        if linkNodeDegree >= linkNodeMaxDegree
            continue;
        end
        
        if Distance < MaxLinkDistance
            fprintf('join cluster[%d,%d],linkIdx = %d\n',row,col,linkIdx);
            amAddCol(linkIdx) = 1;
            amAddRow(linkIdx) = 1;
            
            NewEdgeDelay = 0.5*Distance/100000;
            edgeDelayAddCol(linkIdx) = NewEdgeDelay;
            edgeDelayAddRow(linkIdx) = NewEdgeDelay;
            
            NewLDT = GetLDT(MaxLinkDistance,[nodePos,linkNodePos]);
            ldtAddCol(linkIdx) = NewLDT;
            ldtAddRow(linkIdx) = NewLDT;
            
            vertexStability(nodesNum+1) = NewLDT;
            vertexStability(linkIdx) = vertexStability(linkIdx) + NewLDT;
        end
        break;
    end 
    %�˶δ��������ڴ˴�
    %���½��û�����ӳɹ�������Ӧ��ϢҲ�ܵõ�����
    am = [am,amAddCol];
    am(nodesNum+1,:) = amAddRow;
    
    edgeDelay = [edgeDelay,edgeDelayAddCol];
    edgeDelay(nodesNum+1,:) = edgeDelayAddRow;
    
    ldt = [ldt,ldtAddCol];
    ldt(nodesNum+1,:) = ldtAddRow;
    
    num = nodesNum+1;
    for k = 1:num
        vertexDelay(k) = GetVertexDelay(k,am,edgeDelay);
    end
    
    for k = 1:num
        vertexPriority(k) = GetPriority(k,vertexStability,am,vertexDelay);
    end
    
    AM_{row,col} = am;
    EdgeDelay_{row,col} = edgeDelay;
    VertexDelay_{row,col} = vertexDelay; 
    VertexStability_{row,col} = vertexStability;
    LDT_{row,col} = ldt;
    VertexPriority_{row,col} = vertexPriority;
    IsClusterHead = SetClusterHead(IsClusterHead_,[row,col],vertexPriority);
    
    AM = AM_;
    EdgeDelay = EdgeDelay_;
    VertexDelay = VertexDelay_;
    VertexMaxDegree = VertexMaxDegree_;
    LDT = LDT_;
    VertexStability = VertexStability_;
    VertexPriority = VertexPriority_;
end


