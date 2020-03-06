function [InterClusterInfo] = VertexOutAndIn(ClusterMatrix,MaxLinkDistance,...
    ClusterIdx,nodeIdx,nodePos,AM,VertexMaxDegree,BorderLength,RowBase,ColBase,...
    RowCnt,ColCnt,InterClusterInfo_)
%% 处理的流程为：
%% （1）节点尝试与所在位置的簇内节点相连
%% （2）若成功，则更新节点与相连节点的信息（连接到其它簇，连接到的簇编号，连接到的结点编号）


%%输入参数
%ClusterIdx: 1*2行向量，指明某个簇
%nodeIdx: 数值，指明簇中某个节点
    
    %节点度是否超出
    amOut = AM{ClusterIdx(1),ClusterIdx(2)};
    vertexMaxDegreeOut = VertexMaxDegree{ClusterIdx(1),ClusterIdx(2)};
    infoOut = InterClusterInfo_{ClusterIdx(1),ClusterIdx(2)};
    nodeDegree = sum(amOut(nodeIdx,:)) + size(infoOut{3,nodeIdx},2);
    nodeMaxDegree = vertexMaxDegreeOut(nodeIdx);
    if nodeDegree >= nodeMaxDegree
        InterClusterInfo = InterClusterInfo_;
        return;
    end
    
    %调试
%     fprintf('clusterIdx=[%d,%d],nodeIdx = %d\n',ClusterIdx(1),ClusterIdx(2),nodeIdx);

    row = min(max(ceil(10 * nodePos(3)/(BorderLength*RowBase)),1),RowCnt);
    col = min(max(ceil(10 * nodePos(2)/(BorderLength * ColBase)),1),ColCnt);

    %标记是否已经与所在簇的节点相连,若已相连,
    %(1)若相连的节点与其簇内的大部分节点已脱离，则仍尝试连接簇内的其它结点
    %(2)若相连的节点与其簇内的大部分节点都能连通，则不尝试连接簇内的其它结点了
    if infoOut{1,nodeIdx} == 1
        aLinkClusters = infoOut{2,nodeIdx};
        aLinkIdxs = infoOut{3,nodeIdx};
        linkNum = size(aLinkClusters,2);
        for i = 1:linkNum
            linkClusterIdx = aLinkClusters(:,i);
            if linkClusterIdx(1) == row && linkClusterIdx(2) == col
                aLinkIdx = aLinkIdxs(i);
                [~,MostConnected] = CheckConnected(AM{row,col},aLinkIdx,inf);
                if MostConnected == 1
                    InterClusterInfo = InterClusterInfo_;
                    return;
                end
            end
        end
    end

    amIn = AM{row,col};
    nodesNum = size(amIn,2);
    vertexMaxDegreeIn = VertexMaxDegree{row,col};
    infoIn = InterClusterInfo_{row,col};
    
     %选择距离最近的节点搭建链路(检查链路有效性）
    %（1）若结点linkIdx的度超出最大值，则选择下一结点尝试链接
    % (2）若与结点的距离超出或未超出，则相应地设置参数值，设置完后break
    sortedNodes = SortNodesByDistance(1,nodePos,ClusterMatrix,[row,col],[1:nodesNum]);
    for i = 1:nodesNum  
        linkIdx = sortedNodes(1,i);
        linkNodeDegree = sum(amIn(linkIdx,:))+size(infoIn{3,linkIdx},2);
        linkNodeMaxDegree = vertexMaxDegreeIn(linkIdx);
        Distance = sortedNodes(2,i);
        if linkNodeDegree >= linkNodeMaxDegree
            continue;
        end
        
        if Distance < MaxLinkDistance
%             fprintf('===========================NodeJoin===========================\n');
%             fprintf('join cluster[%d,%d],linkIdx = %d\n',row,col,linkIdx);
            
            infoOut{1,nodeIdx} = 1;
            infoOut{2,nodeIdx} = [infoOut{2,nodeIdx},[row,col]'];
            infoOut{3,nodeIdx} = [infoOut{3,nodeIdx},linkIdx];
            infoIn{1,linkIdx} = 1;
            infoIn{2,linkIdx} = [infoIn{2,linkIdx},[ClusterIdx(1),ClusterIdx(2)]'];
            infoIn{3,linkIdx} = [infoIn{3,linkIdx},nodeIdx];
            InterClusterInfo_{ClusterIdx(1),ClusterIdx(2)} = infoOut;
            InterClusterInfo_{row,col} = infoIn;
            break;
        end
    end
    InterClusterInfo = InterClusterInfo_;
end