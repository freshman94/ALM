%%���ܣ����ص�һ�������룬�����븽��������
function [Idx,IsClusterHead,ClusterMatrix,AM,EdgeDelay,VertexDelay,VertexMaxDegree,...
    LDT,VertexStability,VertexPriority]...
    = ReConnectToOthers(IsClusterHead_,MaxLinkDistance,ClusterMatrix_,ClusterIdx,...
    nodeIdx,AM_,EdgeDelay_,VertexDelay_,VertexMaxDegree_,LDT_,VertexStability_,...
    VertexPriority_,BorderLength,RowBase,ColBase,RowCnt,ColCnt)
 %���ݽڵ�Ĳ�ͬλ�ã���ѡ�������Ĵ��಻ͬ
 %�ڵ��λ�÷�Ϊ�����Ͻǡ����Ͻǡ����½ǡ����½ǣ��ص����ĵ����Ϸ����־�Ϊ���Ͻǣ��Դ����ƣ�
 %���ýڵ����ڴ�Ϊ(i,j),���ѡ�������Ĵ��У�����֤�����ش��ڣ���
 %���Ͻ�(i+1,j)��(i,j-1)����i+1,j-1)
 %���½�(i-1,j)��(i,j-1)��(i-1,j-1)
 %���Ͻ�(i+1,j)��(i,j+1)��(i+1,j+1)
 %���½�(i-1,j)��(i,j+1)��(i-1,j+1)
    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    [isConnected,~,pre_visit] = CheckConnected(am,nodeIdx,inf);
    Idx = nodeIdx+1;
    nodesNum = size(am,2);
    visit = [];
    for i = 1:nodesNum
        if pre_visit(i) == 1
            visit = [visit,i];
        end
    end
    %�ò�������
    if isConnected == 0
       nodesNum = size(am,2);
       partNum = size(visit,2);
       isClusterHead = IsClusterHead_{ClusterIdx(1),ClusterIdx(2)};
       hasClusterHead = 0;
       for i = 1:partNum
           if isClusterHead(visit(i)) == 1
               hasClusterHead = 1;
           end
       end
       vertexMaxDegree = VertexMaxDegree_{ClusterIdx(1),ClusterIdx(2)};
       nodeDegree = sum(am(nodeIdx,:));
       nodeMaxDegree = vertexMaxDegree(nodeIdx);
       %���벿�ֵĴ�С���ܳ����ش�С��һ���Ҳ���������(�����Լ����
       if partNum*2 < nodesNum && hasClusterHead == 0 ...
           && nodeDegree < nodeMaxDegree
           fprintf('===============================PartReConnect==================\n'); 
           Cluster = ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)};
           nodePos = Cluster(:,nodeIdx);
           ReconCluster = GetReconCluster(nodePos,ClusterIdx,BorderLength,RowBase,...
                            ColBase,RowCnt,ColCnt);
           clusterNum = size(ReconCluster,2);
           
           %����
           fprintf('Cluster[%d,%d] nodeIdx = %d part:\n',ClusterIdx(1),ClusterIdx(2),nodeIdx);
           for i = partNum
               fprintf('%d\t',visit(i));
           end
           fprintf('\n');
           
           fprintf('try to ReConnect to Cluster:\n');
           for i = 1:clusterNum
               clusterIdx = ReconCluster{i};
               fprintf('[%d,%d]\t',clusterIdx(1),clusterIdx(2));
           end
           fprintf('\n');
           
           Reconnected = 0;
           for i = 1:clusterNum
               if Reconnected == 1
                   break;
               end
              clusterIdx = ReconCluster{i};
              cluster = ClusterMatrix_{clusterIdx(1),clusterIdx(2)};
              am_ = AM_{clusterIdx(1),clusterIdx(2)};
              vertexMaxDegree_ = VertexMaxDegree_{clusterIdx(1),clusterIdx(2)};
              num = size(cluster,2);
              %���������ӵĴ��еĽڵ㰴����ڵ�nodeIdx�ľ�������
              sortedNodes = SortNodesByDistance(1,nodePos,ClusterMatrix_,...
                            clusterIdx,[1:num]);
              for j = 1:num
                  linkIdx = sortedNodes(1,j);
                  Distance = sortedNodes(2,i);
                  linkNodeDegree = sum(am_(linkIdx,:));
                  linkNodeMaxDegree = vertexMaxDegree_(linkIdx);
                  if linkNodeDegree >= linkNodeMaxDegree
                      continue;
                  end
                  if Distance > MaxLinkDistance
                      break;
                  end
                  %���ӵĲ��֣���Ҫ�����󲿷ֵĽڵ�(50%���ϣ�
                  [~,~,visit_] = CheckConnected(am_,j,inf);
                  if sum(visit_)*2 > num
                      %�ڵ�nodeIdx�������ڵ���ͨ�������룬�������ڴ�
                      fprintf('Reconnect to Cluster[%d,%d],linkIdx = %d\n',...
                          clusterIdx(1),clusterIdx(2),linkIdx);
                      
                      %ͳ���뿪�����б�ű�nodeIdxС��������ȷ��Idxֵ
                      n = 0;
                      for p = 1:partNum
                          if visit(p) < nodeIdx
                              n = n+1;
                          end
                      end
                      Idx = nodeIdx - n;
                      
                      [IsClusterHead_,ClusterMatrix_,AM_,EdgeDelay_,VertexDelay_,...
                        VertexMaxDegree_,LDT_,VertexStability_,VertexPriority_]...
                        = PartOutAndIn(visit,IsClusterHead_,ClusterMatrix_,ClusterIdx,nodeIdx,...
                        clusterIdx,linkIdx,Distance,MaxLinkDistance,AM_,EdgeDelay_,...
                        VertexDelay_,VertexMaxDegree_,LDT_,VertexStability_,...
                        VertexPriority_);
                      Reconnected = 1;
                      break;
                  end                  
              end   
           end
       end
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

function [ReconCluster] = GetReconCluster(nodePos,ClusterIdx,BorderLength,RowBase,...
        ColBase,RowCnt,ColCnt)
     %ȷ���ڵ�����λ�ü���ѡ�������ĸ�����
           %pos��ȡֵ1,2,3,4���ֱ��Ӧ���Ͻǣ����½ǣ����Ͻǣ����½�
           posX = ceil(10 * nodePos(2)/(BorderLength*ColBase)+0.5);
           posY = ceil(10 * nodePos(3)/(BorderLength*RowBase)+0.5);
           row = ClusterIdx(1);
           col = ClusterIdx(2);
           if posX == (row+1) && posY == col
               pos = 1;
               fprintf('���Ͻ�\n');
           elseif posX == row && posY == col
               pos = 2;
               fprintf('���½�\n');
           elseif posX == (row+1) && posY == (col+1)
               pos = 3;
               fprintf('���Ͻ�\n');
           else
               pos = 4;
               fprintf('���½�\n');
           end
           %����
           if pos < 1 || pos > 4
            fprintf('Error!!!! pos = %d\n',pos);
           end
           
           ReconCluster = {};
           switch(pos)
               case 1
                   if row < RowCnt
                       ReconCluster = [ReconCluster,[row+1,col]];
                   end
                   if col > 1
                       ReconCluster = [ReconCluster,[row,col-1]];
                   end
                   if row < RowCnt && col > 1
                       ReconCluster = [ReconCluster,[row+1,col-1]];
                   end
               case 2
                   if row > 1
                       ReconCluster = [ReconCluster,[row-1,col]];
                   end
                   if col > 1
                       ReconCluster = [ReconCluster,[row,col-1]];
                   end
                   if row > 1 && col > 1
                       ReconCluster = [ReconCluster,[row-1,col-1]];
                   end
               case 3
                   if row < RowCnt
                       ReconCluster = [ReconCluster,[row+1,col]];
                   end
                   if col < ColCnt
                       ReconCluster = [ReconCluster,[row,col+1]];
                   end
                   if row < RowCnt && col < ColCnt
                       ReconCluster = [ReconCluster,[row+1,col+1]];
                   end
               case 4
                   if row > 1
                       ReconCluster = [ReconCluster,[row-1,col]];
                   end
                   if col < ColCnt
                       ReconCluster = [ReconCluster,[row,col+1]];
                   end
                   if row > 1 && col < ColCnt
                       ReconCluster = [ReconCluster,[row-1,col+1]];
                   end
           end
end

function [IsClusterHead,ClusterMatrix,AM,EdgeDelay,VertexDelay,...
    VertexMaxDegree,LDT,VertexStability,VertexPriority]...
    = PartOutAndIn(visit,IsClusterHead_,ClusterMatrix_,OutClusterIdx,nodeIdx,...
    InClusterIdx,linkIdx,Distance,MaxLinkDistance,AM_,EdgeDelay_,VertexDelay_,...
    VertexMaxDegree_,LDT_,VertexStability_,VertexPriority_)
    
    Cluster = ClusterMatrix_{OutClusterIdx(1),OutClusterIdx(2)};
    am = AM_{OutClusterIdx(1),OutClusterIdx(2)};
    edgeDelay = EdgeDelay_{OutClusterIdx(1),OutClusterIdx(2)};
    vertexDelay = VertexDelay_{OutClusterIdx(1),OutClusterIdx(2)};
    vertexMaxDegree = VertexMaxDegree_{OutClusterIdx(1),OutClusterIdx(2)};
    ldt = LDT_{OutClusterIdx(1),OutClusterIdx(2)};
    vertexStability = VertexStability_{OutClusterIdx(1),OutClusterIdx(2)};
    vertexPriority = VertexPriority_{OutClusterIdx(1),OutClusterIdx(2)};
    isClusterHead = IsClusterHead_{OutClusterIdx(1),OutClusterIdx(2)};
    
    visitNum = size(visit,2);
    rows = size(Cluster,1);
    part_cluster = zeros(rows,visitNum);
    part_am = zeros(visitNum,visitNum);
    part_edgeDelay = Inf(visitNum,visitNum);
    part_vertexMaxDegree = zeros(1,visitNum);
    part_ldt = zeros(visitNum,visitNum);
    part_vertexStability = zeros(1,visitNum);
    nodePos = Cluster(:,nodeIdx);
    %��ȡ����Ҫ�뿪���ֵ���Ϣ
    for i = 1:visitNum
        if visit(i) == nodeIdx
            newNodeIdx = i;
        end
        part_cluster(:,i) = Cluster(:,visit(i));
        part_vertexMaxDegree(i) = vertexMaxDegree(visit(i));
        part_vertexStability(i) = vertexStability(visit(i));
        for j = i+1:visitNum
            if am(visit(i),visit(j)) == 1
                part_am(i,j) = 1;
                part_am(j,i) = 1;
                part_edgeDelay(i,j) = edgeDelay(visit(i),visit(j));
                part_edgeDelay(j,i) = part_edgeDelay(i,j);
                part_ldt(i,j) = ldt(visit(i),visit(j));
                part_ldt(j,i) = part_ldt(i,j);
            end
        end
    end
    
   %ɾ����Ҫ�뿪���ֵ���Ϣ������Ӻ���ǰɾ��
   visit_sort = sort(visit,'descend');
    for i = 1:visitNum
        am(visit_sort(i),:) = [];
        am(:,visit_sort(i)) = [];

        Cluster(:,visit_sort(i)) = [];
        
        edgeDelay(:,visit_sort(i)) = [];
        edgeDelay(visit_sort(i),:) = [];
        
        ldt(:,visit_sort(i)) = [];
        ldt(visit_sort(i),:) = [];

        vertexDelay(visit_sort(i)) = [];
        vertexMaxDegree(visit_sort(i)) = [];       
        vertexStability(visit_sort(i)) = [];
        vertexPriority(visit_sort(i)) = [];
        
        isClusterHead(visit_sort(i)) = [];
    end
    
    nodesNum = size(am,2);
    for i = 1:nodesNum
        vertexStability(i) = sum(ldt(i,:));
    end
    
    for i = 1:nodesNum
        vertexDelay(i) = GetVertexDelay(i,am,edgeDelay);
    end
    
    for i = 1:nodesNum
        vertexPriority(i) = GetPriority(i,vertexStability,am,vertexDelay);
    end
    
    IsClusterHead = SetClusterHead(IsClusterHead_,OutClusterIdx,vertexPriority);
    
    AM_{OutClusterIdx(1),OutClusterIdx(2)} = am; 
    Cluster(1,:) = [1:nodesNum];
    ClusterMatrix_{OutClusterIdx(1),OutClusterIdx(2)} = Cluster;
    EdgeDelay_{OutClusterIdx(1),OutClusterIdx(2)} = edgeDelay; 
    LDT_{OutClusterIdx(1),OutClusterIdx(2)} = ldt;
    VertexDelay_{OutClusterIdx(1),OutClusterIdx(2)} = vertexDelay;
    VertexMaxDegree_{OutClusterIdx(1),OutClusterIdx(2)} = vertexMaxDegree;
    VertexStability_{OutClusterIdx(1),OutClusterIdx(2)} = vertexStability;
    VertexPriority_{OutClusterIdx(1),OutClusterIdx(2)} = vertexPriority;
    IsClusterHead_{OutClusterIdx(1),OutClusterIdx(2)} = isClusterHead;
    
    %����Ĳ��ּ������ڴ�
    Cluster_ = ClusterMatrix_{InClusterIdx(1),InClusterIdx(2)};
    am_ = AM_{InClusterIdx(1),InClusterIdx(2)};
    edgeDelay_ = EdgeDelay_{InClusterIdx(1),InClusterIdx(2)};
    vertexDelay_ = VertexDelay_{InClusterIdx(1),InClusterIdx(2)};
    vertexMaxDegree_ = VertexMaxDegree_{InClusterIdx(1),InClusterIdx(2)};
    ldt_ = LDT_{InClusterIdx(1),InClusterIdx(2)};
    vertexStability_ = VertexStability_{InClusterIdx(1),InClusterIdx(2)};
    vertexPriority_ = VertexPriority_{InClusterIdx(1),InClusterIdx(2)};  
    inNodesNum = size(am_,2);
    
    for i = 1:visitNum
        part_cluster(1,i) = inNodesNum+i;
        Cluster_(:,inNodesNum+i) = part_cluster(:,i);
    end
    
    amAddRows = zeros(visitNum,inNodesNum);
    amAddCols = zeros(inNodesNum,visitNum);
    amAddCols = [amAddCols;part_am];
    am_ = [am_;amAddRows];
    am_ = [am_,amAddCols];
    am_(linkIdx,newNodeIdx+inNodesNum) = 1;
    am_(newNodeIdx+inNodesNum,linkIdx) = 1;
    
    edgeDelayAddRows = Inf(visitNum,inNodesNum);
    edgeDelayAddCols = Inf(inNodesNum,visitNum);
    edgeDelayAddCols = [edgeDelayAddCols;part_edgeDelay];
    edgeDelay_ = [edgeDelay_;edgeDelayAddRows];
    edgeDelay_ = [edgeDelay_,edgeDelayAddCols];
    NewEdgeDelay = 0.5*Distance/100000;
    edgeDelay_(linkIdx,newNodeIdx+inNodesNum) = NewEdgeDelay;
    edgeDelay_(newNodeIdx+inNodesNum,linkIdx) = NewEdgeDelay;

    vertexMaxDegree_ = [vertexMaxDegree_,part_vertexMaxDegree];
    
    ldtAddRows = zeros(visitNum,inNodesNum);
    ldtAddCols = zeros(inNodesNum,visitNum);
    ldtAddCols = [ldtAddCols;part_ldt];
    ldt_ = [ldt_;ldtAddRows];
    ldt_ = [ldt_,ldtAddCols];
    linkNodePos = Cluster_(:,linkIdx);
    NewLDT = GetLDT(MaxLinkDistance,[nodePos,linkNodePos]);
    ldt_(linkIdx,newNodeIdx+inNodesNum) = NewLDT;
    ldt_(newNodeIdx+inNodesNum,linkIdx) = NewLDT;
    
    vertexStability_ = [vertexStability_,part_vertexStability];
    vertexStability_(linkIdx) = vertexStability_(linkIdx) + NewLDT;
    vertexStability_(newNodeIdx+inNodesNum) = part_vertexStability(newNodeIdx)...
                        +NewLDT;
    num = inNodesNum+visitNum;
    for k = 1:num
        vertexDelay_(k) = GetVertexDelay(k,am_,edgeDelay_);
    end
    
    for k = 1:num
        vertexPriority_(k) = GetPriority(k,vertexStability_,am_,vertexDelay_);
    end
    
    IsClusterHead = SetClusterHead(IsClusterHead_,InClusterIdx,vertexPriority_);
    
    ClusterMatrix_{InClusterIdx(1),InClusterIdx(2)} = Cluster_;
    AM_{InClusterIdx(1),InClusterIdx(2)} = am_;
    EdgeDelay_{InClusterIdx(1),InClusterIdx(2)} = edgeDelay_;
    VertexMaxDegree_{InClusterIdx(1),InClusterIdx(2)} = vertexMaxDegree_;
    LDT_{InClusterIdx(1),InClusterIdx(2)} = ldt_;
    VertexDelay_{InClusterIdx(1),InClusterIdx(2)} = vertexDelay_; 
    VertexStability_{InClusterIdx(1),InClusterIdx(2)} = vertexStability_;
    VertexPriority_{InClusterIdx(1),InClusterIdx(2)} = vertexPriority_;
    
    ClusterMatrix = ClusterMatrix_;
    AM = AM_;
    EdgeDelay = EdgeDelay_;
    VertexDelay = VertexDelay_;
    VertexMaxDegree = VertexMaxDegree_;
    LDT = LDT_;
    VertexStability = VertexStability_;
    VertexPriority = VertexPriority_;    
    
    fprintf('================Out Cluster am=================\n');
    for p = 1:nodesNum
        for q = 1:nodesNum
            fprintf('am(%d,%d) = %d\t',p,q,am(p,q));
        end
        fprintf('\n');
    end
    fprintf('=========================\n');
    
    fprintf('================In Cluster am=================\n');
    for p = 1:num
        for q = 1:num
            fprintf('am_(%d,%d) = %d\t',p,q,am_(p,q));
        end
        fprintf('\n');
    end
    fprintf('=========================\n');
end
