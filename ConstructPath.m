function [edgeTo] = ConstructPath(isClusterHead,am,edgeDelay)
%     fprintf('=================ConstructPath====================\n');
    nodesNum = size(isClusterHead,2);
    for i = 1:nodesNum
        if isClusterHead(i) == 1
            headIdx = i;
            break;
        end
    end
    
    edgeTo = zeros(2,nodesNum);
    distTo = Inf(1,nodesNum);
    distTo(headIdx) = 0;
    IndexMinPQ = [];
    
    IndexMinPQ = PQinsert(IndexMinPQ,headIdx,0);
    PQnum = size(IndexMinPQ,2);
    while PQnum > 0
%         fprintf('IndexMinPQ: ');
%         for i = 1:PQnum
%             rows = size(IndexMinPQ,1);
%             fprintf('[');
%             for j = 1:rows
%                 fprintf('%d ',IndexMinPQ(j,i));
%             end
%             fprintf(']\t');
%         end
%         fprintf('\n');
        [IndexMinPQ,minIdx] = PQdelmin(IndexMinPQ);
%         fprintf('minIdx = %d\n',minIdx);
        [edgeTo,distTo,IndexMinPQ] = relax(am,edgeDelay,minIdx,edgeTo,distTo,...
                                        IndexMinPQ);
        PQnum = size(IndexMinPQ,2);
    end
    
end

function [edgeTo,distTo,IndexMinPQ] = relax(am,edgeDelay,nodeIdx,edgeTo_,distTo_,...
                                        IndexMinPQ_)
    nodesNum = size(am,2);
    for i = 1:nodesNum
        if am(nodeIdx,i) == 1
            if distTo_(i) > distTo_(nodeIdx)+edgeDelay(nodeIdx,i)
                distTo_(i) = distTo_(nodeIdx) + edgeDelay(nodeIdx,i);
                edgeTo_(:,i) = [nodeIdx,i]';
                idx = PQcontain(IndexMinPQ_,i);
                if idx > 0
%                     fprintf('PQ change %d\n',i);
                    IndexMinPQ_ = PQchange(IndexMinPQ_,idx,i,distTo_(i));
                else
%                     fprintf('PQ insert %d\n',i);
                    IndexMinPQ_ = PQinsert(IndexMinPQ_,i,distTo_(i));
                end
            end
        end
    end
    
%     fprintf('edgeTo: ');
%     for i = 1:nodesNum
%         fprintf('%d:[%d->%d]\t',i,edgeTo_(1,i),edgeTo_(2,i));
%     end
%     fprintf('\n');
%     
%     fprintf('distTo: ');
%     for i = 1:nodesNum
%         fprintf('%d:%f\t',i,distTo_(i));
%     end
%     fprintf('\n');
    
    edgeTo = edgeTo_;
    distTo = distTo_;
    IndexMinPQ = IndexMinPQ_;
end

function [IndexMinPQ,minIdx] = PQdelmin(IndexMinPQ_)
    minIdx = IndexMinPQ_(1,1);
    IndexMinPQ_(:,1) = [];
    IndexMinPQ = PQup(IndexMinPQ_);
end

function [IndexMinPQ] = PQinsert(IndexMinPQ_,nodeIdx,priority)
    IndexMinPQ_ = [IndexMinPQ_,[nodeIdx,priority]'];
    IndexMinPQ = PQup(IndexMinPQ_);  
end

function [IndexMinPQ] = PQchange(IndexMinPQ_,idx,nodeIdx,priority)
    IndexMinPQ_(:,idx) = [nodeIdx,priority]';
    IndexMinPQ = PQup(IndexMinPQ_);   
end

function [IndexMinPQ] = PQup(IndexMinPQ_)
    num = size(IndexMinPQ_,2);
    if num > 1
        IndexMinPQ_2 = IndexMinPQ_(2,:);
        IndexMinPQ_2_sort = sort(IndexMinPQ_2);
        pos = find(IndexMinPQ_2==IndexMinPQ_2_sort(1));
        tmp = IndexMinPQ_(:,1);
        IndexMinPQ_(:,1) = IndexMinPQ_(:,pos);
        IndexMinPQ_(:,pos) = tmp;
    end
    
    IndexMinPQ = IndexMinPQ_; 
end

function [idx] = PQcontain(IndexMinPQ,nodeIdx)
    num = size(IndexMinPQ,2);
    idx = 0;
    for i = 1:num
       if IndexMinPQ(1,i) == nodeIdx
           idx = i;
           break;
       end
    end
end