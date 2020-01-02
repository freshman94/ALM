%构建组播树
%随机选择某个簇首为组播源，构建跳数最少的组播树
%(1)src->mediatorA->mediatorB
%loop
%(2)mediatorB->headB
%(3)mediatorB->mediatorC or headB->mediatorC
function [tree,src] = ConstructMTree(ClusterMatrix_,AM_,InterClusterInfo_,...
                    IsClusterHead_,Paths_,EdgeDelay_)
    [RowCnt,ColCnt] = size(IsClusterHead_);
    row = randi(RowCnt,1,1);
    col = randi(ColCnt,1,1);
    num = RowCnt*ColCnt;
    visit = zeros(RowCnt,ColCnt);
    visit(row,col) = 1;
        
    tree = {};
    src = [row,col];
%     src
    %(1)src->mediatorA->mediatorB
    [finish,visit,tree,next_clusters,next_mediators] = ...
        SrcToM(InterClusterInfo_,IsClusterHead_,Paths_,EdgeDelay_,...
        ClusterMatrix_,visit,tree,src,AM_);
    if finish == 1
       return; 
    end
    
    for i = 1:num
        nextNum = size(next_clusters,2);
        next_clusters_ = [];
        next_mediators_ = [];
        for k = 1:nextNum
            next_cluster = next_clusters(:,k);
            next_mediator = next_mediators(k);
            %key
            visit(next_cluster(1),next_cluster(2)) = 1;
            
            %分支（与簇首脱离的部分）
            break_branches = {};
            
            path_ = Paths_{next_cluster(1),next_cluster(2)};
            isClusterHead_ = IsClusterHead_{next_cluster(1),next_cluster(2)};
            nodesNum = size(isClusterHead_,2);
            for p = 1:nodesNum
               if isClusterHead_(p) == 1
                  headIdx_ = p;
                  break;
               end
            end

            %mediator是簇首
            if next_mediator == headIdx_
                tree = [tree;next_cluster.';headIdx_];
                [finish,visit,tree,nextClusters_,nextMediators_]...
                    = SrcToM(InterClusterInfo_,IsClusterHead_,Paths_,EdgeDelay_,...
                    ClusterMatrix_,visit,tree,next_cluster.',AM_);
                if finish == 1
                   break; 
                end
                next_clusters_ = [next_clusters_,nextClusters_];
                next_mediators_ = [next_mediators_,nextMediators_];
                continue;
            end
            %mediator不是簇首
            %(2)mediatorB->headB
            branch_ = next_mediator;
            node = next_mediator;
            %用于指示mediator与Head是否相通
            mToHead = 1;
            for p = 1:nodesNum
                link = path_(:,node);
                up = link(1);
                %mediator与headIdx_不相通，结束
                if up == 0
                    mToHead = 0;
                    break;
                end;
                branch_ = [branch_,up];
                if up == headIdx_
                   tree = [tree;next_cluster.';branch_]; 
                   break;
                end
                node = up;
            end

            if mToHead == 0
               break; 
            end

            %(3)mediatorB->mediatorC or headB->mediatorC
            info_ = InterClusterInfo_{next_cluster(1),next_cluster(2)};
            [row,col] = size(InterClusterInfo_);
            interClusters = zeros(row,col);
            condidates = [];
            branches = {};
            minHop = [];
            minDelay = [];
            am = AM_{next_cluster(1),next_cluster(2)};
            edgeDelay = EdgeDelay_{next_cluster(1),next_cluster(2)};
            isFreshed = 0;
            for j = 1:nodesNum
                if info_{1,j} == 1
                  clusters = info_{2,j};
                  nodes = info_{3,j};
                  cnt = size(clusters,2);
                  for p = 1:cnt
                     cluster = clusters(:,p); 
                     node = nodes(p);

                     isClusterHead_adj = IsClusterHead_{cluster(1),cluster(2)};
                     nodesNum_adj = size(isClusterHead_adj,2);
                     for n = 1:nodesNum_adj
                        if isClusterHead_adj(n) == 1
                           headIdx_adj = n; 
                        end
                     end

                     am_adj = AM_{cluster(1),cluster(2)};
                     [~,~,nodeVisit] = CheckConnected(am_adj,node,inf);
                     %若该节点与其簇首不相通，需将信息亦转发给它，并以该节点为源，
                     %建立其所属的脱离部分的传输路径
                     if nodeVisit(headIdx_adj) == 0
                         %检测是否有重复
                         break_exit = 0;
                         break_branchNum = size(break_branches,1);
                         if break_branchNum ~= 0
                            for n = 1:2:(break_branchNum-1)
                                break_cluster = break_branches{n};
                                break_branch = break_branches{n+1};
                                if size(break_cluster,1) == 1 &&...
                                        break_cluster(1) == cluster(1) && break_cluster(2) == cluster(2)
                                    if ismember(node,break_branch)
                                        break_exit = 1;
                                        break;
                                    end
                                end
                            end
                         end

                         if break_exit == 0
                             break_branches = [break_branches;[next_cluster.';cluster.'];[j,node]];
                             %为脱离部分建立路径
                             if sum(nodeVisit) > 1
                                 edgeDelay_adj = EdgeDelay_{cluster(1),cluster(2)};
                                 nodesNum_adj = size(nodeVisit,2);
                                 for n = 1:nodesNum_adj
                                    if nodeVisit(n) == 1 && n ~= node
                                       minPath = MinPath(node,n,am_adj,edgeDelay_adj);
                                       break_branches = [break_branches;cluster.';minPath];
                                    end
                                 end
                             end
                         end   
                     end

                     if visit(cluster(1),cluster(2)) == 0 && nodeVisit(headIdx_adj) == 1
                         %headB->mediatorC（需要mediator与head能相通)
                         if j == headIdx_ && mToHead == 1
                             interClusters(cluster(1),cluster(2)) = 1;
                             %若已存在与该簇的branch，并且node与condidate(4)可相通
                             condidatesNum = size(condidates,2);
                             for n = 1:condidatesNum
                                condidate = condidates(:,n);
                                if condidate(2) == cluster(1) && condidate(3) == cluster(2)...
                                        && nodeVisit(condidate(4)) == 1
                                    condidates(1,n) = j;
                                    condidates(4,n) = node;
                                    branches{n} = j;
                                    minHop(n,3) = 0;
                                    minDelay(n,3) = 0;
                                    isFreshed = 1;
                                    break;
                                end
                             end
                             if isFreshed == 0 
                                 condidates = [condidates,[j,cluster(1),cluster(2),node]'];
                                 branches = [branches,j];
                                 minHop = [minHop,[cluster(1),cluster(2),0]'];
                                 minDelay = [minDelay,[cluster(1),cluster(2),0]'];
                             end                     
                         else
                             %mediatorB->mediatorC
                             am_adj = AM_{cluster(1),cluster(2)};
                             [~,~,nodeVisit] = CheckConnected(am_adj,node,inf);
                             %若condidates中已存在至该簇的路径
                             %仅当nodeVisit不同时，才添加至condidates
                             canReach = -1;
                             condidatesNum = size(condidates,2);
                             for n = 1:condidatesNum
                                condidate = condidates(:,n);
                                if condidate(2) == cluster(1) && condidate(3) == cluster(2)
                                    canReach = 0;
                                    if nodeVisit(condidate(4)) == 1                                            
                                         canReach = 1;
                                    end
                                end
                             end

                             if canReach == 0
                                 condidates = [condidates,[j,cluster(1),cluster(2),node]'];
                                 branches = [branches,mToIn_path];
                                 minHop = [minHop,[cluster(1),cluster(2),hop]'];
                                 minDelay = [minDelay,[cluster(1),cluster(2),delay]'];
                                continue; 
                             end

                             if interClusters(cluster(1),cluster(2)) == 0
                                 [mToIn_path,hop,delay] = MinPath(next_mediator,j,am,edgeDelay);
                                 if size(mToIn_path,2) == 0
                                    continue;
                                 end
                                 interClusters(cluster(1),cluster(2)) = 1;
                                 condidates = [condidates,[j,cluster(1),cluster(2),node]'];
                                 branches = [branches,mToIn_path];
                                 minHop = [minHop,[cluster(1),cluster(2),hop]'];
                                 minDelay = [minDelay,[cluster(1),cluster(2),delay]'];

                             else
                                 %簇间传输路径有多条，
                                 %(1)选直接与附近簇首相连的(2)选时延最小的
                                 sameInter = 0;
                                 condidatesNum = size(condidates,2);
                                 for n = 1:condidatesNum
                                    condidate = condidates(:,n);
                                    if condidate(1) == j && condidate(2) == cluster(1) &&...
                                            condidate(3) == cluster(2) && condidate(4) ~= headIdx_adj...
                                            && node ~= headIdx_adj
                                        sameInter = 1;
                                        if condidate(4) == headIdx_adj
                                           continue; 
                                        elseif node == headIdx_adj
                                           condidates(4,n) = node; 
                                           continue;
                                        end

                                        cluster_=ClusterMatrix_{next_cluster(1),next_cluster(2)};
                                        interCluster=ClusterMatrix_{cluster(1),cluster(2)};
                                        nodePos = cluster_(:,j);
                                        linkNodePos_1 = interCluster(:,node);
                                        linkNodePos_2 = interCluster(:,condidate(4));
                                        distance_1 = (nodePos(2)-linkNodePos_1(2))^2+...
                                            (nodePos(3)-linkNodePos_1(3))^2;
                                        distance_2 = (nodePos(2)-linkNodePos_2(2))^2+...
                                            (nodePos(3)-linkNodePos_2(3))^2;
                                        if distance_1 < distance_2
                                           condidates(4,n) = node; 
                                           continue;
                                        end
                                    end
                                 end

                                 if sameInter == 1
                                    continue; 
                                 end

                                 [mToIn_path,hop,delay] = MinPath(next_mediator,j,am,edgeDelay);
                                 if size(mToIn_path,2) == 0
                                    continue;
                                 end
                                 minHopNum = size(minHop,2);
                                 for n = 1:minHopNum
                                    if minHop(1,n) == cluster(1) && minHop(2,n) == cluster(2) 
                                        min_idx = n;
                                        break;
                                    end
                                 end
                                 min_hop = minHop(3,min_idx);
                                 if hop < min_hop
                                    condidates(:,min_idx) = [j,cluster(1),cluster(2),node]';
                                    branches{min_idx} = mToIn_path;
                                    minHop(3,min_idx) = hop;
                                    minDelay(3,min_idx) = delay;
                                %跳数相同
                                %(1)选择与附近簇簇首直接相连的（2）选时延最小的
                                 elseif hop == min_hop
                                     min_delay = minDelay(3,min_idx);
                                     if node == headIdx_adj || delay < min_delay
                                        condidates(:,min_idx) = [j,cluster(1),cluster(2),node]';
                                        branches{min_idx} = mToIn_path;
                                        minDelay(3,min_idx) = delay;
                                     end
                                 end
                             end
                         end
                     end
                  end
                end
            end
            condidatesNum = size(condidates,2);
            for n = 1:condidatesNum
                condidate = condidates(:,n);
                branch = branches{n};
                tree = [tree;next_cluster.';branch];
                tree = [tree;[next_cluster.';condidate(2),condidate(3)];[condidate(1),condidate(4)]];
                next_clusters_ = [next_clusters_,[condidate(2),condidate(3)]'];
                next_mediators_ = [next_mediators_,condidate(4)];
                %当condidate(4)与所在簇节点全连通时，所在簇的visit值才置为1
                am_ = AM_{condidate(2),condidate(3)};
                isConnected = CheckConnected(am_,condidate(4),inf);
                if isConnected == 1
                   visit(condidate(2),condidate(3)) = 1; 
                end
            end
            if size(break_branches,1) ~= 0
               tree = [tree;break_branches]; 
            end
        end
        
        if size(next_clusters_,2) == 0
            break;
        end
        next_clusters = next_clusters_;
        next_mediators = next_mediators_;
    end        
end

%src->mediatorA->mediatorB
function [finish,visit,tree,next_clusters,next_mediators]...
        = SrcToM(InterClusterInfo_,IsClusterHead_,Paths_,EdgeDelay_,...
        ClusterMatrix_,visit_,tree_,src,AM_)
    isClusterHead = IsClusterHead_{src(1),src(2)};
    info = InterClusterInfo_{src(1),src(2)};
    path = Paths_{src(1),src(2)};
    edgeDelay = EdgeDelay_{src(1),src(2)};
    nodesNum = size(isClusterHead,2);
    finish = 0;
    for j = 1:nodesNum
        if isClusterHead(j) == 1
            headIdx = j;
        end
    end

    break_branches = {};
    interNodes = [];
    %interNodes是4阶矩阵，每一列包含的信息有:nodeIdx,link_clusterIdx,link_nodeIdx
    for j = 1:nodesNum
       if info{1,j} == 1  
          clusters = info{2,j};
          nodes = info{3,j};
          cnt = size(clusters,2);
          for p = 1:cnt
             cluster = clusters(:,p); 
             node = nodes(p);
             
             isClusterHead_adj = IsClusterHead_{cluster(1),cluster(2)};
             nodesNum_adj = size(isClusterHead_adj,2);
             for n = 1:nodesNum_adj
                if isClusterHead_adj(n) == 1
                   headIdx_adj = n; 
                end
             end
             
             am_adj = AM_{cluster(1),cluster(2)};
             [~,~,nodeVisit] = CheckConnected(am_adj,node,inf);
             %若该节点与其簇首不相通，需将信息亦转发给它，并以该节点为源，
             %建立其所属的脱离部分的传输路径
             if nodeVisit(headIdx_adj)== 0
                 %检测是否有重复
                 break_exit = 0;
                 break_branchNum = size(break_branches,1);
                 if break_branchNum ~= 0
                    for n = 1:2:(break_branchNum-1)
                        break_cluster = break_branches{n};
                        break_branch = break_branches{n+1};
                        if size(break_cluster,1) == 1 &&...
                                break_cluster(1) == cluster(1) && break_cluster(2) == cluster(2)
                            if ismember(node,break_branch)
                                break_exit = 1;
                                break;
                            end
                        end
                    end
                 end
                if break_exit == 0
                     break_branches = [break_branches;[src;cluster.'];[j,node]];
                     %为脱离部分建立路径
                     if sum(nodeVisit) > 1
                         edgeDelay_adj = EdgeDelay_{cluster(1),cluster(2)};
                         nodesNum_adj = size(nodeVisit,2);
                         for n = 1:nodesNum_adj
                            if nodeVisit(n) == 1 && n ~= node
                               minPath = MinPath(node,n,am_adj,edgeDelay_adj);
                               break_branches = [break_branches;cluster.';minPath];
                            end
                         end
                     end
                 end 
             end
             
             if visit_(cluster(1),cluster(2)) == 0 && nodeVisit(headIdx_adj) == 1
                 interNodes = [interNodes,[j,cluster(1),cluster(2),node]']; 
                 continue;
             end
          end
       end
    end
%     interNodes
    interNodesNum = size(interNodes,2);
    if interNodesNum == 0 && size(break_branches,1) == 0
        finish = 1;
        visit = visit_;
        tree = tree_;
        next_clusters = [];
        next_mediators = [];
        return;
    end
    
    %选择与簇首跳数最少的结点（hop同则选时延最小）
    %interClusters表示condidates中是否存在已经连接到某一簇的节点
    [row,col] = size(InterClusterInfo_);
    interClusters = zeros(row,col);
    %minHop是三阶矩阵，表示到达某附近簇的最小跳数
    minHop = [];
    %condidates是四阶矩阵，表示可以与附近簇相连的所有节点
    condidates = [];
    branches = {};
    %统计跳数
    isFreshed = 0;
    for j = 1:interNodesNum
        interNode = interNodes(:,j);
        %簇首与相邻簇相连
        if interNode(1) == headIdx
           interClusters(interNode(2),interNode(3)) = 1; 
           %若已存在与该簇的branch,并且interNode(4)与condidate(4)之间可达，更新
           condidatesNum = size(condidates,2);
           am_adj = AM_{interNode(2),interNode(3)};
           [~,~,nodeVisit] = CheckConnected(am_adj,interNode(4),inf);
           for n = 1:condidatesNum
              condidate = condidates(:,n);
              if condidate(2) == interNode(2) && condidate(3) == interNode(3) &&...
                  nodeVisit(condidate(4)) == 1
                 condidates(1,n) = interNode(1);
                 condidates(4,n) = interNode(4);
                 branches{n} = interNode(1);
                 minHop(3,n) = 0;
                 isFreshed = 1;
                 break;
              end
           end
           
           if isFreshed == 0
               condidates = [condidates,interNode];
               branches = [branches,interNode(1)];
               minHop = [minHop,[interNode(2),interNode(3),0]'];
           end
           continue;
        end
        
        isClusterHead_ = IsClusterHead_{interNode(2),interNode(3)};
        nodesNum_ = size(isClusterHead_,2);
        for n = 1:nodesNum_
           if isClusterHead_(n) == 1
               headIdx_ = n;
           end
        end

        %簇间传输路径有多条，(1)选与附近簇首直接相连的(2)选时延最小的
        sameInter = 0;
        condidatesNum = size(condidates,2);
        for k = 1:condidatesNum
            condidate = condidates(:,k);
            if interNode(1) == condidate(1) && interNode(2) == condidate(2) &&...
                    interNode(3) == condidate(3) && condidate(4) ~= headIdx_...
                    && interNode(4) ~= headIdx_ 
               sameInter = 1;
               if condidate(4) == headIdx_
                   break;
               elseif interNode(4) == headIdx_
                   condidates(4,k) = interNode(4);
                   break;
               end              
               
               cluster = ClusterMatrix_{src(1),src(2)};
               interCluster = ClusterMatrix_{interNode(2),interNode(3)};
               nodePos = cluster(:,interNode(1));
               linkNodePos_1 = interCluster(:,interNode(4));
               linkNodePos_2 = interCluster(:,condidate(4));
               distance_1 = (nodePos(2)-linkNodePos_1(2))^2 +...
                            (nodePos(3)-linkNodePos_1(3))^2;
               distance_2 = (nodePos(2)-linkNodePos_2(2))^2 +...
                            (nodePos(3)-linkNodePos_2(3))^2;
               if distance_1 < distance_2
                   condidates(4,k) = interNode(4);
                   break;
               end
            end
        end
        
        if sameInter == 1
            continue;
        end

        branch = interNode(1); 
        pathNum = size(path,2);
        hop = 0;
        node = interNode(1);
        for p = 1:pathNum
            link = path(:,node);
            up = link(1);
            %节点无法与簇首连通
            if up == 0
               %所有interNode均不能与簇首相通
               if size(condidates,2) == 0 && j == interNodesNum
                   fprintf('==================src broke===============\n');
                   finish  = 1;
                   visit = visit_;
                   tree = tree_;
                   next_clusters = [];
                   next_mediators = [];
                   return;
               end
               break; 
            end
            hop = hop + 1;
            branch = [up,branch];
            if up == headIdx
                %若condidates中已存在至该簇的路径
                %且interNode与condidate(4)不可通时，才添加至condidates
                am_adj = AM_{interNode(2),interNode(3)};
                [~,~,nodeVisit] = CheckConnected(am_adj,interNode(4),inf);
                canReach = -1;
                condidatesNum = size(condidates,2);
                for n = 1:condidatesNum
                   condidate = condidates(:,n);
                   if condidate(2) == interNode(2) && condidate(3) == interNode(3)
                       canReach = 0;
                       if nodeVisit(condidate(4)) == 1
                           canReach = 1;
                       end
                   end
                end
                if canReach == 0
                   condidates = [condidates,interNode];
                   branches = [branches,branch];
                   minHop = [minHop,[interNode(2),interNode(3),hop]'];
                   break; 
                end
                
               if interClusters(interNode(2),interNode(3)) == 0
                   interClusters(interNode(2),interNode(3)) = 1;
                   condidates = [condidates,interNode];
                   branches = [branches,branch];
                   minHop = [minHop,[interNode(2),interNode(3),hop]'];
               else
                   %获取至该附近簇的最小跳数
                   minHopNum = size(minHop,2);
                   for k = 1:minHopNum
                      if minHop(1,k) == interNode(2) && minHop(2,k) == interNode(3)
                         min_hop = minHop(3,k); 
                         hop_idx = k;
                         break;
                      end
                   end
                   %获取相应的condidate和branch
                   condidatesNum = size(condidates,2);
                   for k = 1:condidatesNum
                       condidate = condidates(:,k);
                       if condidate(2)==interNode(2)&&condidate(3)==interNode(3)
                            condidate_idx = k;
                            break;
                       end
                   end
                   if hop < min_hop
                       minHop(3,hop_idx) = hop;
                       %更新condidates,branches
                       condidates(:,condidate_idx) = interNode;
                       branches{condidate_idx} = branch;
                   elseif hop == min_hop
                       %跳数相同，(1)选直接与簇首相连的（2）选时延最小
                       if condidates(4,condidate_idx) == headIdx_
                           break;
                       elseif interNode(4) == headIdx_
                           condidates(:,condidate_idx) = interNode;
                           branches{condidate_idx} = branch;
                           break;
                       end
                       %选时延最小的
                       branchSize = size(branch,2);
                       delay = 0;
                       for k = 1:(branchSize-1)
                          delay = delay+edgeDelay(branch(k),branch(k+1));
                       end
                       con_branch = branches{condidate_idx};
                       con_branchSize = size(con_branch,2);
                       con_delay = 0;
                       for k = 1:(con_branchSize-1)
                           con_delay = con_delay+edgeDelay(con_branch(k),con_branch(k+1));
                       end
                       if con_delay < delay
                          break; 
                       else
                           condidates(:,condidate_idx) = interNode;
                           branches{condidate_idx} = branch;
                           break;
                       end
                   else
                       break;
                   end                       
               end
            end
            node = up;
        end
    end
%         condidates
%         branches
%         branchNum = size(branches,2);
%         for p = 1:branchNum
%            branches{p} 
%         end

%         mediator
    condidatesNum = size(condidates,2);

    next_clusters = [];
    next_mediators = [];
    for k = 1:condidatesNum
       condidate = condidates(:,k);
       branch = branches{k};
       tree_ = [tree_;src;branch];
       next_clusters = [next_clusters,[condidate(2),condidate(3)]'];
       next_mediators = [next_mediators,condidate(4)];
       tree_ = [tree_;[src;[condidate(2),condidate(3)]];[condidate(1),condidate(4)]];
       %当节点condidate(4)与所在簇节点全连通，则所在簇的visit值才置为1
       am_ = AM_{condidate(2),condidate(3)};
       isConnected = CheckConnected(am_,condidate(4),inf);
       if isConnected == 1
           visit_(condidate(2),condidate(3)) = 1;
       end
    end
    
    if size(break_branches,1) ~= 0
       tree_ = [tree_;break_branches];
    end
    
    tree = tree_;
    visit = visit_;
end

%获取src->dst的最短路径，若无路径相通,path为[]
function [path,hop,delay] = MinPath(src,dst,am,edgeDelay)
    [edgeTo,distTo] = Dijkstra(src,am,edgeDelay);
    delay = distTo(dst); 
    nodesNum = size(am,2);
    path = dst;
    hop = 0;
    node = dst;
    for i = 1:nodesNum
        link = edgeTo(:,node);
        hop = hop + 1;
        up = link(1);
        path = [up,path];
        if up == 0
            path = [];
            return;
        end
        if up == src
           return;
        end
        node = up;
    end
end

    
function [edgeTo,distTo] = Dijkstra(node,am,edgeDelay)
    nodesNum = size(am,2); 
    edgeTo = zeros(2,nodesNum);
    distTo = Inf(1,nodesNum);
    distTo(node) = 0;
    IndexMinPQ = [];
    
    IndexMinPQ = PQinsert(IndexMinPQ,node,0);
    PQnum = size(IndexMinPQ,2);
    while PQnum > 0
        [IndexMinPQ,minIdx] = PQdelmin(IndexMinPQ);
        [edgeTo,distTo,IndexMinPQ] = relax(am,edgeDelay,minIdx,edgeTo,distTo,...
                                        IndexMinPQ);
        PQnum = size(IndexMinPQ,2);
    end
    
end

function [edgeTo,distTo,IndexMinPQ] = relax(am,edgeDelay,nodeIdx,edgeTo_,distTo_,...
                                        IndexMinPQ_)
    nodesNum = size(am,2);
    for i = 1:nodesNum
        if am(nodeIdx,i) == 1
            if distTo_(i) > distTo_(nodeIdx)+edgeDelay(nodeIdx,i)
                distTo_(i) = distTo_(nodeIdx) + edgeDelay(nodeIdx,i);
                edgeTo_(:,i) = [nodeIdx,i]';
                idx = PQcontain(IndexMinPQ_,i);
                if idx > 0
                    IndexMinPQ_ = PQchange(IndexMinPQ_,idx,i,distTo_(i));
                else
                    IndexMinPQ_ = PQinsert(IndexMinPQ_,i,distTo_(i));
                end
            end
        end
    end

    edgeTo = edgeTo_;
    distTo = distTo_;
    IndexMinPQ = IndexMinPQ_;
end

function [IndexMinPQ,minIdx] = PQdelmin(IndexMinPQ_)
    minIdx = IndexMinPQ_(1,1);
    IndexMinPQ_(:,1) = [];
    IndexMinPQ = PQup(IndexMinPQ_);
end

function [IndexMinPQ] = PQinsert(IndexMinPQ_,nodeIdx,priority)
    IndexMinPQ_ = [IndexMinPQ_,[nodeIdx,priority]'];
    IndexMinPQ = PQup(IndexMinPQ_);  
end

function [IndexMinPQ] = PQchange(IndexMinPQ_,idx,nodeIdx,priority)
    IndexMinPQ_(:,idx) = [nodeIdx,priority]';
    IndexMinPQ = PQup(IndexMinPQ_);   
end

function [IndexMinPQ] = PQup(IndexMinPQ_)
    num = size(IndexMinPQ_,2);
    if num > 1
        IndexMinPQ_2 = IndexMinPQ_(2,:);
        IndexMinPQ_2_sort = sort(IndexMinPQ_2);
        pos = find(IndexMinPQ_2==IndexMinPQ_2_sort(1));
        tmp = IndexMinPQ_(:,1);
        IndexMinPQ_(:,1) = IndexMinPQ_(:,pos);
        IndexMinPQ_(:,pos) = tmp;
    end
    
    IndexMinPQ = IndexMinPQ_; 
end

function [idx] = PQcontain(IndexMinPQ,nodeIdx)
    num = size(IndexMinPQ,2);
    idx = 0;
    for i = 1:num
       if IndexMinPQ(1,i) == nodeIdx
           idx = i;
           break;
       end
    end
end



