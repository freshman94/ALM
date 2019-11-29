%�����idx����nodeIdx�Ͽ����Ӻ����˵���ͨ��
function [isConnected] = CheckConnected(am,idx,nodeIdx)
%nodeIdxָʾ������idx�Ͽ����ӵĽ��
    nodes = size(am,1);
    visit = zeros(1,nodes);
    fprintf('================================================\n');
    fprintf('DFS visit: ');
    visit = DFS(am,idx,visit,nodes);
    fprintf('\n');
    %����
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
        %�ҵ����ӵĽڵ㣬�����ýڵ�
        if am(idx,i) ~= 0 && visit_(i) == 0 
            visit_(i) = 1;
            fprintf('%d = %d \t',i,visit_(i));
            visit_ = DFS(am,i,visit_,nodes);
        end      
    end
    visit = visit_;
end