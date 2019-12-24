function [IsClusterHead,AM,EdgeDelay,VertexDelay,LDT,VertexStability,VertexPriority]...
    = ReConnect(IsClusterHead_,MaxLinkDistance,ClusterMatrix_,ClusterIdx,nodeIdx,AM_,...
    EdgeDelay_,VertexDelay_,VertexMaxDegree_,LDT_,VertexStability_,VertexPriority_,...
    InterClusterInfo_,RowCnt,ColCnt)
%%若节点脱离，尝试让节点重新连接
%%策略：从am(i,j) == 0的结点中寻找最近的结点，并与其连接
%%(1)若当前结点nodeIdx度超出，break
%%(2)若最近的结点未超出MaxLinkDistance,但度超出，则尝试与下一个最近的结点相连
%%(3)若最近结点已超出MaxLinkDistance，该结点无法重连，break
%%(4)若与最近的结点相连，拓扑通，重连成功，break
%%(5)若与最近结点相连，拓扑仍不通，则使最近的结点重复以上过程
%%(6)若所有结点均已重复以上过程，重连过程结束
    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    [~,MostConnected] = CheckConnected(am,nodeIdx,inf);
    if MostConnected == 0
%         fprintf('=============================ReConnect=========================\n');
        Cluster = ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)};
        edgeDelay = EdgeDelay_{ClusterIdx(1),ClusterIdx(2)};
        ldt = LDT_{ClusterIdx(1),ClusterIdx(2)};
        vertexDelay = VertexDelay_{ClusterIdx(1),ClusterIdx(2)};
        vertexStability = VertexStability_{ClusterIdx(1),ClusterIdx(2)};
        vertexPriority = VertexPriority_{ClusterIdx(1),ClusterIdx(2)};
        nodesNum = size(am,1);
        vertexMaxDegree = VertexMaxDegree_{ClusterIdx(1),ClusterIdx(2)};
        info = InterClusterInfo_{ClusterIdx(1),ClusterIdx(2)};
%         for m = 1:nodesNum
%             for n = 1:nodesNum
%                 fprintf('am(%d,%d) = %d\t',m,n,am(m,n));
%             end
%             fprintf('\n');
%         end
%         fprintf('\n');
        for i = 1:nodesNum
            nodeDegree = sum(am(nodeIdx,:));
            nodeMaxDegree = vertexMaxDegree(nodeIdx);
            if nodeDegree >= nodeMaxDegree
                break;
            end
            NotConnectedNodes = [];
            for k = 1:nodesNum
                if am(nodeIdx,k) == 0
                    NotConnectedNodes = [NotConnectedNodes,k];
                end
            end            
            
            nodePos = Cluster(:,nodeIdx);
            sortedNodes = SortNodesByDistance(0,nodePos,ClusterMatrix_,ClusterIdx,...
                NotConnectedNodes);
            Distance = inf;
            
            isRemoved = 0;
            %表明是否可以有结点可连(仅当距离超出时，才不可连）
            isAccess = 0;
            NotConnectedNodesNum = size(sortedNodes,2);
            for p = 1:NotConnectedNodesNum
                linkIdx = sortedNodes(1,p);
                Distance = sortedNodes(2,p);
                if Distance > MaxLinkDistance
                    break;
                end
                linkNodeDegree = sum(am(linkIdx,:));
                linkNodeMaxDegree = vertexMaxDegree(linkIdx);
                if linkNodeDegree >= linkNodeMaxDegree
                    [~,~,visit_pre ] = CheckConnected(am,linkIdx,inf);
                    %尝试断开linkIdx的某条无影响的簇内链路
                    for k = 1:nodesNum
                       if am(linkIdx,k) == 1
                          am(linkIdx,k) = 0;
                          am(k,linkIdx) = 0;
                          [~,~,visit_cur] = CheckConnected(am,linkIdx,inf);
                          %找到某条无影响的链路
                          if isequal(visit_pre,visit_cur) 
%                               fprintf('Cluster[%d,%d]\n',ClusterIdx(1),ClusterIdx(2));
%                               fprintf('Reconnect:Remove inner Link(%d,%d)\n',linkIdx,k);
%                               fprintf('nodeIdx %d Reconnect to linkIdx %d\n',nodeIdx,linkIdx);
                              isRemoved = 1;
                              isAccess = 1;
                              break; 
                          else
                              am(linkIdx,k) = 1;
                              am(k,linkIdx) = 1;
                          end
                       end
                    end
                    if isRemoved == 1
                        break;
                    end
                else
                    isRemoved = 1;
                    isAccess = 1;
                    break;
                end
            end
            %尝试断开某条无影响的簇间链路（不管能否找到无影响，一定会断一条）
            if isRemoved == 0 && isAccess == 1
                for p = 1:NotConnectedNodesNum
                    linkIdx = sortedNodes(1,p);
                    linkNodeDegree = sum(am(linkIdx,:));
                    linkNodeMaxDegree = vertexMaxDegree(linkIdx);
                    if linkNodeDegree >= linkNodeMaxDegree
                        [~,cluster_visit_pre] = ClusterConnected(InterClusterInfo_,...
                            RowCnt,ColCnt);
                        if info{1,linkIdx} == 1
                            linkClusters = info{2,linkIdx};
                            linkIdxs = info{3,linkIdx};
                            linkNum = size(linkClusters,2);
                            for k = 1:linkNum
                                linkCluster = linkClusters(:,k);
                                Cluster_linkIdx = linkIdxs(k);
                                linkClusters(:,k) = [];
                                linkIdxs(k) = [];
                                info{2,linkIdx} = linkClusters;
                                info{3,linkIdx} = linkIdxs;
                                InterClusterInfo_{ClusterIdx(1),ClusterIdx(2)} = info;
                                
                                info_ = InterClusterInfo_{linkCluster(1),linkCluster(2)};
                                linkClusters_ = info_{2,Cluster_linkIdx};
                                linkIdxs_ = info_{3,Cluster_linkIdx};
                                linkNum_ = size(linkClusters_,2);
                                for m = 1:linkNum_
                                   if isequal(linkClusters_(:,m),ClusterIdx.') &&...
                                           linkIdxs_(m) == linkIdx
                                      linkClusters_(:,m) = [];
                                      linkIdxs_(m) = [];
                                      info_{2,Cluster_linkIdx} = linkClusters_;
                                      info_{3,Cluster_linkIdx} = linkIdxs_;
                                      InterClusterInfo_{linkCluster(1),linkCluster(2)} = info_;
                                   end
                                end
                                
                                [~,cluster_visit_cur] = ClusterConnected(InterClusterInfo_,...
                                    RowCnt,ColCnt);
                                if isequal(cluster_visit_pre,cluster_visit_cur) || (k == linkNum)
%                                     for m = 1:nodesNum
%                                         for n = 1:nodesNum
%                                             fprintf('am(%d,%d) = %d\t',m,n,am(m,n));
%                                         end
%                                         fprintf('\n');
%                                     end
%                                     fprintf('\n');
%                                     fprintf('Reconnect:Remove Inter Link([%d,%d]%d,[%d,%d]%d)\n',...
%                                         ClusterIdx(1),ClusterIdx(2),linkIdx,linkCluster(1),...
%                                         linkCluster(2),Cluster_linkIdx);
%                                     fprintf('nodeIdx %d Reconnect to linkIdx %d\n',nodeIdx,linkIdx);
                                    if size(linkClusters,2) == 0
                                        info{1,linkIdx} = [];
                                        InterClusterInfo_{ClusterIdx(1),ClusterIdx(2)} = info;
                                    end
                                    
                                    if size(linkClusters_,2) == 0
                                        info_{1,Cluster_linkIdx} = [];
                                        InterClusterInfo_{linkCluster(1),linkCluster(2)} = info_;
                                    end
                                    
                                    isRemoved = 1;
                                    break;
                                else
                                    info{2,linkIdx} = [linkCluster,info{2,linkIdx}];
                                    info{3,linkIdx} = [Cluster_linkIdx,info{3,linkIdx}];
                                    InterClusterInfo_{ClusterIdx(1),ClusterIdx(2)} = info;
                                    info_{2,Cluster_linkIdx} = [ClusterIdx.',info_{2,Cluster_linkIdx}];
                                    info_{3,Cluster_linkIdx} = [linkIdx,info_{3,Cluster_linkIdx}];
                                    InterClusterInfo_{linkCluster(1),linkCluster(2)} = info_;
                                end  
                            end
                        end
                    end
                end
            end
            
            %添加新链路
            if isRemoved == 1 && isAccess == 1
                am(nodeIdx,linkIdx) = 1;
                am(linkIdx,nodeIdx) = 1;
%                 fprintf('nodeIdx = %d,linkIdx = %d\n',nodeIdx,linkIdx);
%                 for m = 1:nodesNum
%                     for n = 1:nodesNum
%                         fprintf('am(%d,%d) = %d\t',m,n,am(m,n));
%                     end
%                     fprintf('\n');
%                 end
%                 fprintf('\n');
                NewEdgeDelay = 0.5*Distance/100000;
                edgeDelay(nodeIdx,linkIdx) = NewEdgeDelay;
                edgeDelay(linkIdx,nodeIdx) = edgeDelay(nodeIdx,linkIdx);

                for k = 1:nodesNum
                    vertexDelay(k) = GetVertexDelay(k,am,edgeDelay);
                end

                ldt(nodeIdx,linkIdx) = GetLDT(MaxLinkDistance,Cluster(:,[nodeIdx,linkIdx]));
                ldt(linkIdx,nodeIdx) = ldt(nodeIdx,linkIdx);

                vertexStability(nodeIdx) = sum(ldt(nodeIdx,:));
                vertexStability(linkIdx) = sum(ldt(linkIdx,:));

                for k = 1:nodesNum
                    vertexPriority(i) = GetPriority(i,vertexStability,...
                        am,vertexDelay,info);
                end
                IsClusterHead_ = SetClusterHead(IsClusterHead_,ClusterIdx,vertexPriority);

                AM_{ClusterIdx(1),ClusterIdx(2)} = am;
                EdgeDelay_{ClusterIdx(1),ClusterIdx(2)} = edgeDelay;
                LDT_{ClusterIdx(1),ClusterIdx(2)} = ldt;
                VertexDelay_{ClusterIdx(1),ClusterIdx(2)} = vertexDelay;
                VertexStability_{ClusterIdx(1),ClusterIdx(2)} = vertexStability;
                VertexPriority_{ClusterIdx(1),ClusterIdx(2)} = vertexPriority;

                isConnected = CheckConnected(am,nodeIdx,inf);
                if isConnected == 1
                    break;
                else
                    nodeIdx = linkIdx;
                end
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