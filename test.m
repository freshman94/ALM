
clc;
clf;
[ClusterMatrix]=...
    cluster_topology(1e2,40,1e5,5e1,[3,5],[0,3],[0,2*pi]);
hold off;
    

% a = cell(3,1);
% a{1} = [1,1];
% a{2} = [2,2,2];
% a{3} = [3,3,3,3];
% 
% a(2) = [];
% a
