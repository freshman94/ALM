function [paths,visit,stack] = GetPaths(src,des,am,paths_,visit_,stack_)
%pathsΪ��Ԫ����
%visit_Ϊ1*N��������NΪ���н������
    visit_(src) = 1;
    %��ջ
    stack_ = [stack_,src];
    nodesNum = size(visit_,2);
    for i = 1:nodesNum
        if src == des
            %����
%             fprintf('==================Getpaths=================\n');
%             fprintf('src = %d,des = %d\n',stack_(1),des);
%             stackSize = size(stack_,2);
%             for j = 1:stackSize;
%                 fprintf('%d -> ',stack_(j));
%             end
%             fprintf('\n=======================================\n');
            
            %�ҵ�·�����洢����
            paths_ = [paths_,stack_];
            %��ջ
            stackSize = size(stack_,2);
            stack_ = stack_(1:stackSize-1);
            visit_(src) = 0;
            break;
        end
        if am(src,i) == 1 && visit_(i) == 0
            [paths_,visit_,stack_] = GetPaths(i,des,am,paths_,visit_,stack_);
        end
        if i == nodesNum
            stackSize = size(stack_,2);
            stack_ = stack_(1:stackSize-1);
            visit_(src) = 0;
        end
    end
    paths = paths_;
    visit = visit_;
    stack = stack_;
end