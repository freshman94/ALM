%%���ܣ������븽��������
function [IsClusterHead,VertexPriority,InterClusterInfo]...
    = ReConnectToOthers(IsClusterHead_,MaxLinkDistance,ClusterMatrix_,ClusterIdx,...
    nodeIdx,AM_,VertexDelay_,VertexMaxDegree_,VertexStability_,VertexPriority_,BorderLength,RowBase,...
    ColBase,RowCnt,ColCnt,InterClusterInfo_)
 %���ݽڵ�Ĳ�ͬλ�ã���ѡ�������Ĵ��಻ͬ
 %�ڵ��λ�÷�Ϊ�����Ͻǡ����Ͻǡ����½ǡ����½ǣ��ص����ĵ����Ϸ����־�Ϊ���Ͻǣ��Դ����ƣ�
 %���ýڵ����ڴ�Ϊ(i,j),���ѡ�������Ĵ��У�����֤�����ش��ڣ���
 %���Ͻ�(i+1,j)��(i,j-1)����i+1,j-1)
 %���½�(i-1,j)��(i,j-1)��(i-1,j-1)
 %���Ͻ�(i+1,j)��(i,j+1)��(i+1,j+1)
 %���½�(i-1,j)��(i,j+1)��(i-1,j+1)
   am = AM_{ClusterIdx(1),ClusterIdx(2)};
   vertexMaxDegree = VertexMaxDegree_{ClusterIdx(1),ClusterIdx(2)};
   info = InterClusterInfo_{ClusterIdx(1),ClusterIdx(2)};
   nodeDegree = sum(am(nodeIdx,:))+size(info{3,nodeIdx},2);
   nodeMaxDegree = vertexMaxDegree(nodeIdx);
   if nodeDegree < nodeMaxDegree
%        fprintf('===============================ConnectToAdj==================\n'); 
       Cluster = ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)};
       nodePos = Cluster(:,nodeIdx);
       ReconCluster = GetReconCluster(nodePos,ClusterIdx,BorderLength,RowBase,...
                        ColBase,RowCnt,ColCnt);
       clusterNum = size(ReconCluster,2);

       %����
%        fprintf('cluster[%d,%d] nodeIdx = %d\n',ClusterIdx(1),ClusterIdx(2),nodeIdx);
%        fprintf('try to Connect to Cluster:\n');
%        for i = 1:clusterNum
%            clusterIdx = ReconCluster{i};
%            fprintf('[%d,%d]\t',clusterIdx(1),clusterIdx(2));
%        end
%        fprintf('\n');

       for i = 1:clusterNum
          clusterIdx = ReconCluster{i};
          cluster = ClusterMatrix_{clusterIdx(1),clusterIdx(2)};
          am_ = AM_{clusterIdx(1),clusterIdx(2)};
          vertexMaxDegree_ = VertexMaxDegree_{clusterIdx(1),clusterIdx(2)};
          linkInfo_ = InterClusterInfo_{clusterIdx(1),clusterIdx(2)};
          num = size(cluster,2);
          %���������ӵĴ��еĽڵ㰴����ڵ�nodeIdx�ľ�������
          sortedNodes = SortNodesByDistance(1,nodePos,ClusterMatrix_,...
                        clusterIdx,[1:num]);
          isRemoved = 0;
          for j = 1:num
              linkIdx = sortedNodes(1,j);
              Distance = sortedNodes(2,j);
              if Distance > MaxLinkDistance
                  break;
              end
              linkNodeDegree = sum(am_(linkIdx,:))+size(linkInfo_{3,linkIdx},2);
              linkNodeMaxDegree = vertexMaxDegree_(linkIdx);
              [~,~,visit_pre] = CheckConnected(am_,linkIdx,inf);
              %�����������Ľڵ�p���Ѵ����ֵ��p���ԶϿ�ĳ����Ӱ�����·
              if linkNodeDegree >= linkNodeMaxDegree
                  %���ԶϿ���Ӱ��ĳ��������·
                  for k = 1:num
                    if am_(linkIdx,k) == 1
                        am_(linkIdx,k) = 0;
                        am_(k,linkIdx) = 0;
                        [~,~,visit_cur] = CheckConnected(am_,linkIdx,inf);
                        %�ҵ�ĳ����Ӱ����·
                        if isequal(visit_pre,visit_cur)
%                             fprintf('ReconnectToOthers:Remove inner Link(%d,%d)\n',...
%                                 linkIdx,k);
                            isRemoved = 1;
                            break;
                        end
                    end
                  end
                  if isRemoved == 1
                      break;
                  %���ԶϿ�ĳ����Ӱ��ؼ���·
                  else
                    if linkInfo_{1,linkIdx} == 1
                        [~,visit_cluster_pre] = ClusterConnected(InterClusterInfo_,...
                            RowCnt,ColCnt);
                        linkClusters = linkInfo_{2,linkIdx};
                        linkIdxs = linkInfo_{3,linkIdx};
%                         linkInfo_
%                         linkIdx
%                         clusterIdx
%                         linkClusters
%                         linkIdxs
                        linkNum = size(linkClusters,2);
                        for k = 1:linkNum
                            linkCluster = linkClusters(:,k);
                            cluster_linkIdx = linkIdxs(k);
                            linkClusters(:,k) = [];
                            linkIdxs(k) = [];
                            linkInfo_{2,linkIdx} = linkClusters;
                            linkInfo_{3,linkIdx} = linkIdxs;
                            InterClusterInfo_{clusterIdx(1),clusterIdx(2)} = linkInfo_;
                             
                            info_ = InterClusterInfo_{linkCluster(1),linkCluster(2)};
                            linkClusters_ = info_{2,cluster_linkIdx};
                            linkIdxs_ = info_{3,cluster_linkIdx};
%                             info_
%                             cluster_linkIdx
%                             linkClusters_
%                             linkIdxs_
                            linkNum_ = size(linkClusters_,2);
                            for m = 1:linkNum_
                                if isequal(linkClusters_(:,m),clusterIdx.') &&...
                                       linkIdxs_(m) == linkIdx
%                                   fprintf('success\n');
                                  linkClusters_(:,m) = [];
                                  linkIdxs_(m) = [];
                                  info_{2,cluster_linkIdx} = linkClusters_;
                                  info_{3,cluster_linkIdx} = linkIdxs_;
                                  InterClusterInfo_{linkCluster(1),linkCluster(2)} = info_;
                                  break;
                                end
                            end

                             [~,visit_cluster_cur] = ClusterConnected(InterClusterInfo_,...
                                 RowCnt,ColCnt);
                             if isequal(visit_cluster_pre,visit_cluster_cur)
%                                  fprintf('ReconnectToOthers:Remove Inter Link([%d,%d]%d,[%d,%d]%d)\n',...
%                                         clusterIdx(1),clusterIdx(2),linkIdx,...
%                                         linkCluster(1),linkCluster(2),cluster_linkIdx);
                                 if size(linkInfo_{2,linkIdx},2) == 0
                                   linkInfo_{1,linkIdx} = 0;
                                   InterClusterInfo_{clusterIdx(1),clusterIdx(2)} = linkInfo_;
                                 end
                                 
                                 if size(linkClusters_,2) == 0
                                    info_{1,cluster_linkIdx} = 0;
                                    InterClusterInfo_{linkCluster(1),linkCluster(2)} = info_;
                                  end
                                 
                                  break;
                             else
                                 linkInfo_{2,linkIdx} = [linkCluster,linkInfo_{2,linkIdx}];
                                 linkInfo_{3,linkIdx} = [cluster_linkIdx,linkInfo_{3,linkIdx}];
                                 InterClusterInfo_{clusterIdx(1),clusterIdx(2)} = linkInfo_;
                                 info_{2,cluster_linkIdx} = [clusterIdx.',info_{2,cluster_linkIdx}];
                                 info_{3,cluster_linkIdx} = [linkIdx,info_{3,cluster_linkIdx}];
                                 InterClusterInfo_{linkCluster(1),linkCluster(2)} = info_;
                             end
%                              linkClusters
                        end
                    end
                  end
              end

              %���ӵĲ��֣���Ҫ�����󲿷ֵĽڵ�(50%���ϣ�
              if sum(visit_pre)*2 > num
%                    fprintf('[%d,%d]%d Connect to [%d,%d]%d\n',ClusterIdx(1),...
%                        ClusterIdx(2),nodeIdx,clusterIdx(1),clusterIdx(2),linkIdx);
                  [IsClusterHead_,VertexPriority_,InterClusterInfo_]...
                    = ConnectToAdj(IsClusterHead_,ClusterIdx,nodeIdx,clusterIdx,...
                    linkIdx,AM_,VertexDelay_,VertexStability_,VertexPriority_,...
                    InterClusterInfo_);
              end                  
          end
          nodeDegree = sum(am(nodeIdx,:))+size(info{3,nodeIdx},2);
          if nodeDegree >= nodeMaxDegree
              break;
          end
       end
   end

    VertexPriority = VertexPriority_;
    IsClusterHead = IsClusterHead_;
    InterClusterInfo = InterClusterInfo_;
end

function [ReconCluster] = GetReconCluster(nodePos,ClusterIdx,BorderLength,RowBase,...
        ColBase,RowCnt,ColCnt)
     %ȷ���ڵ�����λ�ü���ѡ�������ĸ�����
           %pos��ȡֵ1,2,3,4���ֱ��Ӧ���Ͻǣ����½ǣ����Ͻǣ����½�
           posRow = ceil(10 * nodePos(3)/(BorderLength*RowBase)+0.5);
           posCol = ceil(10 * nodePos(2)/(BorderLength*ColBase)+0.5);
           row = ClusterIdx(1);
           col = ClusterIdx(2);
           if posRow == (row+1) && posCol == col
               pos = 1;
%                fprintf('���Ͻ�\n');
           elseif posRow == row && posCol == col
               pos = 2;
%                fprintf('���½�\n');
           elseif posRow == (row+1) && posCol == (col+1)
               pos = 3;
%                fprintf('���Ͻ�\n');
           else
               pos = 4;
%                fprintf('���½�\n');
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

function [IsClusterHead,VertexPriority,InterClusterInfo]...
    = ConnectToAdj(IsClusterHead_,ClusterIdx,nodeIdx,AdjClusterIdx,linkIdx,...
    AM_,VertexDelay_,VertexStability_,VertexPriority_,InterClusterInfo_)
    
    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    vertexDelay = VertexDelay_{ClusterIdx(1),ClusterIdx(2)};
    vertexStability = VertexStability_{ClusterIdx(1),ClusterIdx(2)};
    vertexPriority = VertexPriority_{ClusterIdx(1),ClusterIdx(2)};
    info = InterClusterInfo_{ClusterIdx(1),ClusterIdx(2)};
    
    %��ǽ���Ƿ�������ڴصĽ������������������
    %(1)�������Ľڵ�������ڵĴ󲿷ֽڵ������룬���Գ������Ӵ��ڵ��������
    %(2)�������Ľڵ�������ڵĴ󲿷ֽڵ㶼����ͨ���򲻳������Ӵ��ڵ����������
    if info{1,nodeIdx} == 1
        aLinkClusters = info{2,nodeIdx};
        aLinkIdxs = info{3,nodeIdx};
        linkNum = size(aLinkClusters,2);
        for i = 1:linkNum
            linkClusterIdx = aLinkClusters(:,i);
            if isequal(linkClusterIdx,AdjClusterIdx.')
                aLinkIdx = aLinkIdxs(i);
                [~,MostConnected] = CheckConnected(AM_{AdjClusterIdx(1),...
                    AdjClusterIdx(2)},aLinkIdx,inf);
                if aLinkIdx == linkIdx || MostConnected == 1
                    IsClusterHead = IsClusterHead_;
                    VertexPriority = VertexPriority_;
                    InterClusterInfo = InterClusterInfo_;
                    return;
                end
            end
        end
    end
    
%     fprintf('[%d,%d]%d Connect to [%d,%d] %d\n',ClusterIdx(1),ClusterIdx(2),nodeIdx,...
%               AdjClusterIdx(1),AdjClusterIdx(2),linkIdx);
    
    info{1,nodeIdx} = 1;
    info{2,nodeIdx} = [info{2,nodeIdx},[AdjClusterIdx(1),AdjClusterIdx(2)]'];
    info{3,nodeIdx} = [info{3,nodeIdx},linkIdx];
%     info
%     info{2,nodeIdx}
%     info{3,nodeIdx}
    
    nodesNum = size(am,2);
    for i = 1:nodesNum
        vertexPriority(i) = GetPriority(i,vertexStability,am,vertexDelay,info);
    end
    
    IsClusterHead_ = SetClusterHead(IsClusterHead_,ClusterIdx,vertexPriority);
    VertexPriority_{ClusterIdx(1),ClusterIdx(2)} = vertexPriority;
    InterClusterInfo_{ClusterIdx(1),ClusterIdx(2)} = info;
    
    %�������ڴص�info
    am_ = AM_{AdjClusterIdx(1),AdjClusterIdx(2)};
    vertexDelay_ = VertexDelay_{AdjClusterIdx(1),AdjClusterIdx(2)};
    vertexStability_ = VertexStability_{AdjClusterIdx(1),AdjClusterIdx(2)};
    vertexPriority_ = VertexPriority_{AdjClusterIdx(1),AdjClusterIdx(2)};  
    linkInfo_ = InterClusterInfo_{AdjClusterIdx(1),AdjClusterIdx(2)};  
    num = size(am_,2);
    
    linkInfo_{1,linkIdx} = 1;
    linkInfo_{2,linkIdx} = [linkInfo_{2,linkIdx},[ClusterIdx(1),ClusterIdx(2)]'];
    linkInfo_{3,linkIdx} = [linkInfo_{3,linkIdx},nodeIdx];
%     linkInfo_
%     linkInfo_{2,linkIdx}
%     linkInfo_{3,linkIdx}
    
    for k = 1:num
        vertexPriority_(k) = GetPriority(k,vertexStability_,am_,vertexDelay_,linkInfo_);
    end
    
    IsClusterHead = SetClusterHead(IsClusterHead_,AdjClusterIdx,vertexPriority_);
    VertexPriority_{AdjClusterIdx(1),AdjClusterIdx(2)} = vertexPriority_;
    InterClusterInfo_{AdjClusterIdx(1),AdjClusterIdx(2)} = linkInfo_;
    
    VertexPriority = VertexPriority_;   
    InterClusterInfo = InterClusterInfo_;
end
