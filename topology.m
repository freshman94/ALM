%���°�Ľ�Salam����������������㷨ͨ��MATLABԴ��
%{
������Ϊ���°�Դ�룬Դ����ɾ�����ܻ��Ư������������ͼƬ���㷨�Ľ�˵�����£�
1.ʹ��K��ֵ������ƽڵ�ֲ������ܣ�ʹ�ò���������������ͨ�Ժ;����Ը���
2.�����������������ݷḻ����������·�ķ��á�ʱ�ӡ������ڵ�ķ��á�ʱ�ӡ�ʱ�Ӷ�����������
3.��·ʱ�ӵ��ڽڵ�����������֮�����٣����ӷ���ʵ�����
%}

function [Sxy,AM,EdgeCost,EdgeDelay,EdgeBandWide,VertexCost,VertexDelay,VertexDelayJitter,VertexPacketLoss]=topology(BorderLength,NodeAmount, ...
Alpha,Beta,PlotIf,EdgeCostDUB,EdgeBandWideDUB,VertexCostDUB,VertexDelayDUB,VertexDelayJitterDUB,VertexPacketLossDUB)
%% ��������б�
%BorderLength������������������ı߳�����λ��km
%NodeAmount������������ڵ�ĸ���
%Alpha����������������������AlphaԽ�󣬶̱���Գ��ߵı���Խ��
%Beta����������������������BetaԽ�󣬱ߵ��ܶ�Խ��
%PlotIf���������Ƿ���������ͼ�����Ϊ1����ͼ�����򲻻�ͼ
%EdgeCostDUB����������·���õĿ��Ʋ�����1*2���洢��·���õ��½���Ͻ�
%EdgeBandWideDUB����������·����Ŀ��Ʋ�����1*2���洢�½���Ͻ�
%VertexCostDUB���������ڵ���õĿ��Ʋ�����1*2,�洢�ڵ���õ��½���Ͻ�
%VertexDelayDUB���������ڵ�ʱ�ӵĿ��Ʋ�����1*2���ڴ��ڵ�ʱ�ӵ��½���Ͻ�
%VertexDelayJitterDUB���������ڵ�ʱ�Ӷ����Ŀ��Ʋ�����1*2���洢�ڵ�ʱ�Ӷ������½���Ͻ�
%VertexPacketLossDUB���������ڵ㶪���ʵĿ��Ʋ�����1*2,�洢�ڵ㶪���ʵ��½�

%%  �������
%Sxy��������3*N�ľ��󣬸��зֱ����ڴ洢�ڵ����ţ������꣬������ľ���
%AM��������0 1�洢����AM(i,j)=1��ʾ������i��j������ߣ�N*N
%EdgeCost����������·���þ���N*N
%EdgeDelay����������·ʱ�Ӿ���N*N
%EdgeBandWide����������·�������N*N 
%VertexCost���������ڵ��������,1*N
%VertexDelay���������ڵ�ʱ��������1*N
%VertexDelayJitter���������ڵ�ʱ�Ӷ�������,1*N
%VertexPacketLoss���������ڵ㶪��������,1*N
%%�Ƽ�������������� 
%BorderLength=1000;NodeAmount=25;Alpha=100000000;Beta=200000000000;
%PlotIf=1;EdgeCostDUB=[2,5];EdgeBandWideDUB=[30,1000];VertexCostDUB=[2,4];
%VertexDelayDUB=1e-4*[5,20];VertexDelayJitterDUB=1e-4*[3,8];
%VertexPacketLossDUB=1e-4*[0,500]
%%
%������ʼ��
NN = 10*NodeAmount;
SSxy = zeros(NN,2);
%���������������������ѡȡNN���ڵ�
for i = 1:NN
    SSxy(i,1) = BorderLength*rand;
    SSxy(i,2) = BorderLength*rand;
end

[IDX,C] = kmeans(SSxy,NodeAmount);Sxy = [[1:NodeAmount]',C]';
%����������С�����˳������Ϊÿһ���ڵ���
temp = Sxy;
Sxy2 = Sxy(2,:);    %�ڶ���Ԫ�أ����ڵ�����꣩
Sxy2_sort = sort(Sxy2);
for i = 1:NodeAmount
    pos = find(Sxy2==Sxy2_sort(i)); %����Sxy2������������Ԫ�ص�λ��
    if length(pos)>1
        error('������ϣ������ԣ�');
    end
    temp(1,i) = i;
    temp(2,i) = Sxy(2,pos);
    temp(3,i) = Sxy(3,pos);
end
Sxy = temp;
%���������ʼ��
AM = zeros(NodeAmount,NodeAmount);
EdgeCost = zeros(NodeAmount,NodeAmount);
EdgeDelay = zeros(NodeAmount,NodeAmount);
EdgeBandWide = zeros(NodeAmount,NodeAmount);
VertexCost  = zeros(1,NodeAmount);
VertexDelay = zeros(1,NodeAmount);
VertexDelayJitter = zeros(1,NodeAmount);
VertexPacketLoss  = zeros(1,NodeAmount);
for i = 1:(NodeAmount-1)
    for j = (i+1):NodeAmount
        Distance =( (Sxy(2,i)-Sxy(2,j))^2+(Sxy(3,i)-Sxy(3,j))^2)^0.5;
        P = Beta*exp(-Distance^5/(Alpha*BorderLength));
        if P>rand
            AM(i,j) = 1;
            AM(j,i) = 1;
            EdgeDelay(i,j) = 0.5*Distance/100000;
            EdgeDelay(j,i) = EdgeDelay(i,j);
            EdgeCost(i,j) = EdgeCostDUB(1)+(EdgeCostDUB(2)-EdgeCostDUB(1))*rand;
            EdgeCost(j,i)=EdgeCost(i,j);
            EdgeBandWide(i,j) = EdgeBandWideDUB(1)+(EdgeBandWideDUB(2)-EdgeBandWideDUB(1))*rand;
            EdgeBandWide(j,i)=EdgeBandWide(i,j);
        else
            EdgeDelay(i,j) = inf;
            EdgeDelay(j,i) = inf;
            EdgeCost(i,j) = inf;
            EdgeCost(j,i) = inf;
            EdgeBandWide(i,j) = inf;
            EdgeBandWide(j,i) = inf;
        end
    end
end
for i = 1:NodeAmount
    VertexCost(i) = VertexCostDUB(1)+(VertexCostDUB(2)-VertexCostDUB(1))*rand;
    VertexDelay(i) = VertexDelayDUB(1)+(VertexDelayDUB(2)-VertexDelayDUB(1))*rand;
    VertexDelayJitter(i) = VertexDelayJitterDUB(1)+(VertexDelayJitterDUB(2)-VertexDelayJitterDUB(1))*rand;
    VertexPacketLoss(i) = VertexPacketLossDUB(1)+(VertexPacketLossDUB(2)-VertexPacketLossDUB(1))*rand;
end
Net_plot(BorderLength,NodeAmount,Sxy,EdgeCost,PlotIf);
end

%���ڻ����������˵ĺ���
function Net_plot(BorderLength,NodeAmount,Sxy,EdgeCost,PlotIf)
%���ڵ�
if PlotIf == 1
    plot(Sxy(2,:),Sxy(3,:),'ko','MarkerEdgeColor','b','MarkerFaceColor','g','MarkerSize',5);    %'ko'����ɫԲȦ��'MarkerEdgeColor'����ǵı߿���ɫ��'MarkerFaceColor'����ǵ���ɫ
    %����ͼ����ʾ��Χ
    xlim([0,BorderLength]);
    ylim([0,BorderLength]);
    hold on;    %��ʾ����ԭͼ���޸�
    %�ڵ�����
    for i = 1:NodeAmount
        Str = int2str(i);
        text(Sxy(2,i)+BorderLength/100,Sxy(3,i)+BorderLength/100,Str,'FontName','Times New Roman','FontSize',12);
        hold on;
    end
end
%����
if PlotIf == 1
    for i = 1:(NodeAmount-1)
        for j = (i+1):NodeAmount
            if isinf(EdgeCost(i,j)) == 0
                plot([Sxy(2,i),Sxy(2,j)],[Sxy(3,i),Sxy(3,j)]);
                hold on;
            end
        end
    end
end
if PlotIf == 1
    xlabel('x (km)','FontName','Times New Roman','FontSize',12);
    ylabel('y (km)','FontName','Times New Roman','FontSize',12);
end
end