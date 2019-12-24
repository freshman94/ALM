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

%�Ӵ�[1,1]��ʼ��Ѱ�Ҹôؿ��������Ĵ�p���Ҵ�p��δ������
%�Ӵ�p��ʼ���ظ����Ϲ��̣�ֱ�����дر�����

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