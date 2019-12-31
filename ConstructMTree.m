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
            %mediator不是簇首
            %(2)mediatorB->headB
            else
                branch_ = next_mediator;
                node = next_mediator;
                %用于指示mediator与Head是否相通
                mToHead = 1;
                for p = 1:nodesNum
                    link = path_(:,node);
                    up = link(1);
                    %mediator与headIdx_不相通
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
                         
                         if visit(cluster(1),cluster(2)) == 0
                             %headB->mediatorC（需要mediator与head能相通)
                             if j == headIdx_ && mToHead == 1
                                 %若已存在与该簇的branch，并且node与condidate(4)可相通
                                 am_adj = AM_{cluster(1),cluster(2)};
                                 [isConnected,~,nodeVisit] = CheckConnected(am_adj,node,inf);
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
                                 
                                 if isConnected == 1
                                    interClusters(cluster(1),cluster(2)) = 1;
                                 end
                                 fprintf('headB->mediatorC:isConnected=%d\n',isConnected);
                             else
                                 %mediatorB->mediatorC
                                 if interClusters(cluster(1),cluster(2)) == 0
                                     [mToIn_path,hop,delay] = MinPath(next_mediator,j,am,edgeDelay);
                                     if size(mToIn_path,2) == 0
                                        continue;
                                     end
                                     am_adj = AM_{cluster(1),cluster(2)};
                                     [isConnected,~,nodeVisit] = CheckConnected(am_adj,node,inf);
                                     if isConnected == 1
                                        interClusters(cluster(1),cluster(2)) = 1;
                                     end
                                     %若condidates中已存在至该簇的路径
                                     %仅当nodeVisit不同时，才添加至condidates
                                     isAdded = 0;
                                     condidatesNum = size(condidates,2);
                                     for n = 1:condidatesNum
                                        condidate = condidates(:,n);
                                        if condidate(2) == cluster(1) && condidate(3) == cluster(2)
                                            if nodeVisit(condidate(4)) == 0
                                                 condidates = [condidates,[j,cluster(1),cluster(2),node]'];
                                                 branches = [branches,mToIn_path];
                                                 minHop = [minHop,[cluster(1),cluster(2),hop]'];
                                                 minDelay = [minDelay,[cluster(1),cluster(2),delay]'];
                                            end
                                            isAdded = 1;
                                            break;
                                        end
                                     end
                                     
                                     if isAdded == 0
                                         condidates = [condidates,[j,cluster(1),cluster(2),node]'];
                                         branches = [branches,mToIn_path];
                                         minHop = [minHop,[cluster(1),cluster(2),hop]'];
                                         minDelay = [minDelay,[cluster(1),cluster(2),delay]'];
                                     end
                                     
                                 else
                                     %簇间传输路径有多条，
                                     %(1)选直接与附近簇首相连的(2)选时延最小的
                                     condidatesNum = size(condidates,2);
                                     for n = 1:condidatesNum
                                        condidate = condidates(:,n);
                                        if condidate(1) == j && condidate(2) == cluster(1) &&...
                                                condidate(3) == cluster(2) && condidate(4) ~= headIdx_adj...
                                                && node ~= headIdx_adj
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
                                               fprintf('==========shorter interLink==========\n'); 
                                               condidates(4,n) = node; 
                                            else
                                                continue;
                                            end
                                        end
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
                    visit(condidate(2),condidate(3)) = 1;
                    tree = [tree;next_cluster.';branch];
                    tree = [tree;[next_cluster.';condidate(2),condidate(3)];[condidate(1),condidate(4)]];
                    next_clusters_ = [next_clusters_,[condidate(2),condidate(3)]'];
                    next_mediators_ = [next_mediators_,condidate(4)];
                    %若condidate(4)不能与所在簇节点全连通，则所在簇的visit值应置为0
                    am_ = AM_{condidate(2),condidate(3)};
                    isConnected = CheckConnected(am_,condidate(4),inf);
                    if isConnected == 0
                       visit(condidate(2),condidate(3)) = 0; 
                    end
                end              
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
             if visit_(cluster(1),cluster(2)) == 0
                 interNodes = [interNodes,[j,cluster(1),cluster(2),node]']; 
                 continue;
             end
          end
       end
    end
%     interNodes
    interNodesNum = size(interNodes,2);
    if interNodesNum == 0
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
           %若已存在与该簇的branch,并且interNode(4)与condidate(4)之间可达，更新
           condidatesNum = size(condidates,2);
           am_adj = AM_{interNode(2),interNode(3)};
           [isConnected,~,nodeVisit] = CheckConnected(am_adj,interNode(4),inf);
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
           
           %即将与其相连的附近簇只有是全连通的，才可以被标记为已访问
           if isConnected == 1
              interClusters(interNode(2),interNode(3)) = 1;  
           end
           fprintf('SrcToM head->interNode:isConnected=%d\n',isConnected);
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
        condidatesNum = size(condidates,2);
        for k = 1:condidatesNum
            condidate = condidates(:,k);
            if interNode(1) == condidate(1) && interNode(2) == condidate(2) &&...
                    interNode(3) == condidate(3) && condidate(4) ~= headIdx_...
                    && interNode(4) ~= headIdx_                
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
                   fprintf('=========shorter interLink==============\n');
               end
            end
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
               if interClusters(interNode(2),interNode(3)) == 0
                   %若condidates中已存在至该簇的路径
                   %仅当nodeVisit不同时，才添加至condidates
                   am_adj = AM_{interNode(2),interNode(3)};
                   [isConnected,~,nodeVisit] = CheckConnected(am_adj,interNode(4),inf);
                   if isConnected == 1
                        interClusters(interNode(2),interNode(3)) = 1;
                   end
                   isAdded = 0;
                   condidatesNum = size(condidates,2);
                   for n = 1:condidatesNum
                       condidate = condidates(:,n);
                       if condidate(2) == interNode(2) && condidate(3) == interNode(3)
                           if nodeVisit(condidate(4)) == 0
                               condidates = [condidates,interNode];
                               branches = [branches,branch];
                               minHop = [minHop,[interNode(2),interNode(3),hop]'];
                           end
                           isAdded = 1;
                           break;
                       end
                   end
                   
                   if isAdded == 0
                       condidates = [condidates,interNode];
                       branches = [branches,branch];
                       minHop = [minHop,[interNode(2),interNode(3),hop]'];
                   end
                   
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



