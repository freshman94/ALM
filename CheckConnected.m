%�����idx����nodeIdx�Ͽ����Ӻ����˵���ͨ��
%am(idx,nodeIdx) == 0
%����������idx����ͨ�ԣ���nodeIdx��ֵΪinf����
function [isConnected,MostConnected,visit] = CheckConnected(am,idx,nodeIdx)
%%���룺nodeIdxָʾ������idx�Ͽ����ӵĽ��
%%�����isConnected��ʾ�Ƿ�ȫ��ͨ
%%MostConnected��ʾ�Ƿ�󲿷���ͨ
%%visitΪ���idx���Է��ʵĽ�㼯��
    nodesNum = size(am,1);
    visit = zeros(1,nodesNum);
%     fprintf('================================================\n');
%     fprintf('DFS visit: ');
    visit = DFS(am,idx,visit,nodesNum);
%     fprintf('\n');
    %����
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
        %�ҵ����ӵĽڵ㣬�����ýڵ�
        if am(idx,i) ~= 0 && visit_(i) == 0 
            visit_(i) = 1;
%             fprintf('%d = %d \t',i,visit_(i));
            visit_ = DFS(am,i,visit_,nodesNum);
        end      
    end
    visit = visit_;
end