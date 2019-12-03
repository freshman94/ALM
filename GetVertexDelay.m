function [VertexDelay] = GetVertexDelay(nodeIdx,linkIdx,ClusterIdx,AM_,VertexDelay_,newEdgeDelay)
%%AM_需要是还未修改的
%%适用于同一个簇内的节点重连
    fprintf('===================GetVertexDelay=====================\n');
    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    nodeDegree = sum(am(nodeIdx,:));
    linkNodeDegree = sum(am(linkIdx,:));
    
    vertexDelay = VertexDelay_{ClusterIdx(1),ClusterIdx(2)};
    nodeDelayPre = vertexDelay(nodeIdx);
    linkNodeDelayPre = vertexDelay(linkIdx);
    
    if isnan(nodeDelayPre) || isinf(nodeDelayPre)
        nodeDelayPre = 0;
    end
    if isnan(linkNodeDelayPre) || isinf(linkNodeDelayPre)
        linkNodeDelayPre = 0;
    end
    
    nodeDelayCur = (nodeDelayPre*nodeDegree + newEdgeDelay)/(nodeDegree+1);
    linkNodeDelayCur = (linkNodeDelayPre*linkNodeDegree+newEdgeDelay)/(linkNodeDegree+1);
    
    fprintf('nodeIdx = %d,nodeDegree = %d,nodeDelayPre = %f,nodeDelayCur = %f\n',...
        nodeIdx,nodeDegree,nodeDelayPre,nodeDelayCur);
    fprintf('linkIdx = %d,linkNodeDegree = %d,linkNodeDelayPre = %f,nodeDelayCur = %f\n',...
        linkIdx,linkNodeDegree,linkNodeDelayPre,linkNodeDelayCur);
    
    vertexDelay(nodeIdx) = nodeDelayCur;
    vertexDelay(linkIdx) = linkNodeDelayCur;
    
    VertexDelay_{ClusterIdx(1),ClusterIdx(2)} = vertexDelay;
    VertexDelay = VertexDelay_;
end