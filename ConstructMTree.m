%构建组播树
%随机选择某个簇首为组播源，构建跳数最少的组播树
%(1)src->mediatorA->mediatorB
%loop
%(2)mediatorB->headB
%(3)mediatorB->mediatorC or headB->mediatorC
function [tree,src] = ConstructMTree(AM_,InterClusterInfo_,IsClusterHead_,Paths_,EdgeDelay_)
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
    [finish,visit,tree,next_cluster,mediator] = SrcToM(InterClusterInfo_,IsClusterHead_,...
            Paths_,EdgeDelay_,visit,tree,src);
    if finish == 1
       return; 
    end
    
    for i = 1:num
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
        if mediator == headIdx_
            tree = [tree;next_cluster.';headIdx_];
            [finish,visit,tree,next_cluster,mediator] = SrcToM(InterClusterInfo_,...
                IsClusterHead_,Paths_,EdgeDelay_,visit,tree,next_cluster.');
            if finish == 1
               break; 
            end   
        %mediator不是簇首
        %(2)mediatorB->headB
        else
            branch_ = mediator;
            node = mediator;
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
            interNodes = [];
            headInter = 0;
            for j = 1:nodesNum
                if (headInter == 0) && (info_{1,j} == 1)
                  clusters = info_{2,j};
                  nodes = info_{3,j};
                  cnt = size(clusters,2);
                  for p = 1:cnt
                     cluster = clusters(:,p); 
                     node = nodes(p);
                     if visit(cluster(1),cluster(2)) == 0
                         %headB->mediatorC（需要mediator与head能相通)
                         if j == headIdx_ && mToHead == 1
                             tree = [tree;[next_cluster.';cluster.'];[headIdx_,node]];
                             next_cluster = cluster;
                             mediator = node;
                             headInter = 1;
                             break;
                         end
                         interNodes = [interNodes,j]; 
                         break;
                     end
                  end
                end
            end
            %mediatorB->mediatorC
            if headInter == 0
                interNodesNum = size(interNodes,2);
                branch = [];
                minHop = inf;
                minDelay = inf;
                minNode = 0;
                for p = 1:interNodesNum
                    interNode = interNodes(p);
                    am = AM_{next_cluster(1),next_cluster(2)};
                    edgeDelay = EdgeDelay_{next_cluster(1),next_cluster(2)};
                    [mToIn_path,hop,delay] = MinPath(mediator,interNode,am,edgeDelay);
                    if size(mToIn_path,2) == 0
                        continue;
                    else
                        if hop < minHop
                           minHop = hop;
                           minDelay = delay;
                           branch = mToIn_path;
                           minNode = interNode;
                        elseif hop == minHop
                            if delay < minDelay
                                minDelay = delay;
                                branch = mToIn_path;
                                minNode = interNode;
                            end
                        end
                    end   
                end
                %无路径相通
                if size(branch,2) == 0
                   return; 
                end
                tree = [tree;next_cluster.';branch];
                
                clusters = info_{2,minNode};
                nodes = info_{3,minNode};
                cnt = size(clusters,2);
                for p = 1:cnt
                   cluster = clusters(:,p);
                   node = nodes(p);
                   if visit(cluster(1),cluster(2)) == 0
                      visit(cluster(1),cluster(2)) = 1;
                      tree = [tree;[next_cluster.';cluster.'];[minNode,node]];
                      next_cluster = cluster;
                      mediator = node;
                      break;
                   end
                   %簇间无路径相通
                   if p == cnt
                       return;
                   end
                end
            end
        end
    end
             
end

%src->mediatorA->mediatorB
function [finish,visit,tree,next_cluster,next_mediator] = SrcToM(InterClusterInfo_,IsClusterHead_,...
                        Paths_,EdgeDelay_,visit_,tree_,src)
    tree_ = [tree_;src];
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
    for j = 1:nodesNum
       if info{1,j} == 1  
          clusters = info{2,j};
          cnt = size(clusters,2);
          for p = 1:cnt
             cluster = clusters(:,p); 
             if visit_(cluster(1),cluster(2)) == 0
                 interNodes = [interNodes,j]; 
                 break;
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
        next_cluster = [];
        next_mediator = 0;
        return;
    else 
        %选择与簇首跳数最少的结点（hop同则选时延最小）
        minHop = inf;
        condidates = [];
        branches = {};
        %统计跳数
        for j = 1:interNodesNum
            interNode = interNodes(j);
            %簇首与相邻簇相连
            if interNode == headIdx
               condidates = interNode;
               branches = {interNode};
               break;
            end
            branch = interNode; 
            pathNum = size(path,2);
            hop = 0;
            node = interNode;
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
                       next_cluster = [];
                       next_mediator = 0;
                       return;
                   end
                   break; 
                end
                hop = hop + 1;
                branch = [up,branch];
                if up == headIdx
                   if hop < minHop  
                      minHop = hop;
                      condidates = [];
                      branches = {};
                      condidates(1) = interNode;
                      branches{1} = branch;

                    elseif hop == minHop
                      condidates = [condidates,interNode]; 
                      branches = [branches,branch];
                   end
                   break; 
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
        %跳数相同,(1)选直接与簇首相连的（2）选时延最小
        condidatesNum = size(condidates,2);
        mediator = 0;
        if condidatesNum > 1
            minDelay = inf;
            seq = 0;
            for j = 1:condidatesNum
                condidate = condidates(j);
                branch = branches{j};
                %判断是否与相邻簇的簇首直接相连
                clusters = info{2,condidate};
                nodes = info{3,condidate};
                cnt = size(clusters,2);
                for m = 1:cnt
                    cluster = clusters(:,m);
                    if visit_(cluster(1),cluster(2)) == 0
                        node = nodes(m);
                        isClusterHead_ = IsClusterHead_{cluster(1),cluster(2)};
                        nodesNum_ = size(isClusterHead_,2);
                        for n = 1:nodesNum_
                           if isClusterHead_(n) == 1
                               headIdx_ =  n;
                               break;
                           end
                        end
                        if node == headIdx_
                            tree_ = [tree_;branches{j}];
                            next_cluster = cluster;
                            next_mediator = node;
                            tree_ = [tree_;[src;cluster.'];[condidate,node]];
                            visit_(cluster(1),cluster(2)) = 1;
                            visit = visit_;
                            tree = tree_;
                            fprintf('==============Link Head=====================\n');
                            return;
                        end
                    end
                end
                    
                branchSize = size(branch,2);
                delay = 0;
                for p = 1:(branchSize-1)
                    delay = delay+edgeDelay(branch(p),branch(p+1));
                end
                if delay < minDelay
                   minDelay = delay;
                   mediator = condidate;
                   seq = j;
                end
            end
            tree_ = [tree_;branches{seq}];
        else
           mediator = condidates(1); 
           tree_ = [tree_;branches{1}];
        end
%         mediator
        clusters = info{2,mediator};
        nodes = info{3,mediator};
        cnt = size(clusters,2);
        for p = 1:cnt
            cluster = clusters(:,p);
           if visit_(cluster(1),cluster(2)) == 0
              visit_(cluster(1),cluster(2)) = 1;
              node = nodes(p);
              next_cluster = cluster;
              next_mediator = node;
              tree_ = [tree_;[src;cluster.'];[mediator,node]];
              break;
           end
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



