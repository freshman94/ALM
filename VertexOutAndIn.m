function [InterClusterInfo] = VertexOutAndIn(ClusterMatrix,MaxLinkDistance,...
    ClusterIdx,nodeIdx,nodePos,AM,VertexMaxDegree,BorderLength,RowBase,ColBase,...
    RowCnt,ColCnt,InterClusterInfo_)
%% ���������Ϊ��
%% ��1���ڵ㳢��������λ�õĴ��ڽڵ�����
%% ��2�����ɹ�������½ڵ��������ڵ����Ϣ�����ӵ������أ����ӵ��Ĵر�ţ����ӵ��Ľ���ţ�


%%�������
%ClusterIdx: 1*2��������ָ��ĳ����
%nodeIdx: ��ֵ��ָ������ĳ���ڵ�
    
    %�ڵ���Ƿ񳬳�
    amOut = AM{ClusterIdx(1),ClusterIdx(2)};
    vertexMaxDegreeOut = VertexMaxDegree{ClusterIdx(1),ClusterIdx(2)};
    infoOut = InterClusterInfo_{ClusterIdx(1),ClusterIdx(2)};
    nodeDegree = sum(amOut(nodeIdx,:)) + size(infoOut{3,nodeIdx},2);
    nodeMaxDegree = vertexMaxDegreeOut(nodeIdx);
    if nodeDegree >= nodeMaxDegree
        InterClusterInfo = InterClusterInfo_;
        return;
    end
    
    %����
%     fprintf('clusterIdx=[%d,%d],nodeIdx = %d\n',ClusterIdx(1),ClusterIdx(2),nodeIdx);

    row = min(max(ceil(10 * nodePos(3)/(BorderLength*RowBase)),1),RowCnt);
    col = min(max(ceil(10 * nodePos(2)/(BorderLength * ColBase)),1),ColCnt);

    %����Ƿ��Ѿ������ڴصĽڵ�����,��������,
    %(1)�������Ľڵ�������ڵĴ󲿷ֽڵ������룬���Գ������Ӵ��ڵ��������
    %(2)�������Ľڵ�������ڵĴ󲿷ֽڵ㶼����ͨ���򲻳������Ӵ��ڵ����������
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
    
     %ѡ���������Ľڵ���·(�����·��Ч�ԣ�
    %��1�������linkIdx�Ķȳ������ֵ����ѡ����һ��㳢������
    % (2��������ľ��볬����δ����������Ӧ�����ò���ֵ���������break
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