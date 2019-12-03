%������·������ʱ��
function [LDT] = GetLDT(MaxLinkDistance,linkNodePos)
    fprintf('=========================GetLDT=========================\n');
    b = linkNodePos(2,1) - linkNodePos(2,2);
    d = linkNodePos(3,1) - linkNodePos(3,2);
    a = linkNodePos(4,1)*cos(linkNodePos(5,1)) - linkNodePos(4,2)*cos(linkNodePos(5,2));
    c = linkNodePos(4,1)*sin(linkNodePos(5,1)) - linkNodePos(4,2)*sin(linkNodePos(5,2));
    
    LDT = (-1*(a*b+c*d)+sqrt( (a.^2+ c.^2)*MaxLinkDistance.^2 - (a*d-b*c).^2))...
            / (a.^2+ c.^2);
    fprintf('speed = [%f,%f],direc = [%f,%f]\n',linkNodePos(4,1),linkNodePos(4,2),...
        linkNodePos(5,1),linkNodePos(5,2));
    %��������ڵ���ٶȺͷ���һ�£�LDTΪinf����Ӱ������ֵ�ļ��㣬���תΪһ���ܴ����������
    if isinf(LDT) || isnan(LDT)
        LDT = 10e10;
    end
end