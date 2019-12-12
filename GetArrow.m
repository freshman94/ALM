function [x,y] = GetArrow(src,des,BorderLength)
    alpha = atan2(des(2)-src(2),des(1)-src(1));
    fprintf('alpha = %f\n',alpha);
    direc_1 = alpha + pi/6;
    direc_2 = alpha - pi/6;
    
    mid = [(src(1)+des(1))/2,(src(2)+des(2))/2];
    fprintf('mid:[%d,%d]\n',mid(1),mid(2));
    beg_1 = [mid(1)-BorderLength/100*cos(direc_1),mid(2)-BorderLength/100*sin(direc_1)];
    beg_2 = [mid(1)-BorderLength/100*cos(direc_2),mid(2)-BorderLength/100*sin(direc_2)];
    fprintf('beg_1:[%d,%d]\n',beg_1(1),beg_1(2));
    fprintf('beg_2:[%d,%d]\n',beg_2(1),beg_2(2));
    x = [beg_1(1),mid(1),beg_2(1)];
    y = [beg_1(2),mid(2),beg_2(2)]; 
end