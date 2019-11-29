%检测结点idx与结点nodeIdx断开链接后，拓扑的连通性
function [isConnected] = CheckConnected(am,idx,nodeIdx)
%nodeIdx指示即将与idx断开链接的结点
    nodes = size(am,1);
    visit = zeros(1,nodes);
    fprintf('================================================\n');
    fprintf('DFS visit: ');
    visit = DFS(am,idx,visit,nodes);
    fprintf('\n');
    %调试
    fprintf('visit: ');
    for i = 1:nodes
        fprintf('%d = %d\t',i,visit(i));
    end
    fprintf('\n================================================\n');
    
    for i = 1:nodes
        if  i ~= nodeIdx && visit(i) == 0
            isConnected = 0;
            return;
        end
    end
    isConnected = 1;
end

function [visit] = DFS (am,idx,visit_,nodes)
    visit_(idx) = 1;
    
    for i = 1:nodes
        %找到连接的节点，遍历该节点
        if am(idx,i) ~= 0 && visit_(i) == 0 
            visit_(i) = 1;
            fprintf('%d = %d \t',i,visit_(i));
            visit_ = DFS(am,i,visit_,nodes);
        end      
    end
    visit = visit_;
end