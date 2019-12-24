%检测结点idx与结点nodeIdx断开链接后，拓扑的连通性
%am(idx,nodeIdx) == 0
%若仅是想检测idx的连通性，则将nodeIdx赋值为inf即可
function [isConnected,MostConnected,visit] = CheckConnected(am,idx,nodeIdx)
%%输入：nodeIdx指示即将与idx断开链接的结点
%%输出：isConnected表示是否全连通
%%MostConnected表示是否大部分连通
%%visit为结点idx可以访问的结点集合
    nodesNum = size(am,1);
    visit = zeros(1,nodesNum);
%     fprintf('================================================\n');
%     fprintf('DFS visit: ');
    visit = DFS(am,idx,visit,nodesNum);
%     fprintf('\n');
    %调试
%     fprintf('visit: ');
%     for i = 1:nodesNum
%         fprintf('%d = %d\t',i,visit(i));
%     end
%     fprintf('\n================================================\n');
    include = 0;
    if nodeIdx < nodesNum
        include = 1;
    end
    isConnected = 0;
    MostConnected = 0;
    visitNum = sum(visit);
    if nodesNum == (visitNum - include)
        isConnected = 1;
    end
    MostPercent = 0.7;
    if visitNum >= floor(MostPercent*nodesNum) 
        MostConnected = 1;
    end
end

function [visit] = DFS (am,idx,visit_,nodesNum)
    visit_(idx) = 1;
    
    for i = 1:nodesNum
        %找到连接的节点，遍历该节点
        if am(idx,i) ~= 0 && visit_(i) == 0 
            visit_(i) = 1;
%             fprintf('%d = %d \t',i,visit_(i));
            visit_ = DFS(am,i,visit_,nodesNum);
        end      
    end
    visit = visit_;
end