function [isConnected,visit] = ClusterConnected(InterClusterInfo,RowCnt,ColCnt)
    visit = zeros(RowCnt,ColCnt);
    visit = DFS(InterClusterInfo,[1,1],visit);
    for i = 1:RowCnt
        for j = 1:ColCnt
            if visit(i,j) == 0
                isConnected = 0;
                return;
            end
        end
    end
    isConnected = 1;
end

%从簇[1,1]开始，寻找该簇可以相连的簇p，且簇p还未被访问
%从簇p开始，重复以上过程，直至所有簇被访问

function [visit] = DFS(InterClusterInfo,ClusterIdx,visit_)
    info = InterClusterInfo{ClusterIdx(1),ClusterIdx(2)};
    visit_(ClusterIdx(1),ClusterIdx(2)) = 1;
    nodesNum = size(info,2);
    for i = 1:nodesNum
        if info{1,i} == 1
            linkClusters = info{2,i};
            linkNum = size(linkClusters,2);
            for j = 1:linkNum  
                linkClusterIdx = linkClusters(:,j);
                if visit_(linkClusterIdx(1),linkClusterIdx(2)) == 0
                    visit_(linkClusterIdx(1),linkClusterIdx(2)) = 1;
                    visit_ = DFS(InterClusterInfo,linkClusterIdx,visit_);
                end
            end
        end
    end
    visit = visit_;
end