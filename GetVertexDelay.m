function [delay] = GetVertexDelay(nodeIdx,am,edgeDelay)
%更新节点的时延，am,edgeDelay需要是更新后的

    fprintf('=====================GetVertexDelay==================\n');
    nodesNum = size(am,2);
    node2nodeDelay = zeros(1,nodesNum);
    for i = 1:nodesNum
        if i ~= nodeIdx
            paths = {};
            visit = zeros(1,nodesNum);
            stack = [];
            [paths,~,~] = GetPaths(nodeIdx,i,am,paths,visit,stack);
            pathNum = size(paths,2);
            if pathNum == 0
                node2nodeDelay(i) = 10e10;
            else
                pathDelay = 0;
                for m = 1:pathNum
                    path = paths{m};
                    pathNodesNum = size(path,2);
                    for n = 1:pathNodesNum-1
                        pathDelay = pathDelay + edgeDelay(path(n),path(n+1));
                    end
%                     fprintf('pathDelay(%d,%d) = %f\t',nodeIdx,i,pathDelay);
                end
%                 fprintf('\n');
                node2nodeDelay(i) = pathDelay/pathNum;
%                 fprintf('node2nodeDelay(%d,%d) = %f\t',nodeIdx,i,node2nodeDelay(i));
            end
        end
    end
    delay = sum(node2nodeDelay)/(nodesNum-1);                    
end
