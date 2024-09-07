% find truncation point
function [leftInx, rightInx] = ftrup(rirBright)
    lenInx = 1:length(rirBright);

    [pks1, locs1] = findpeaks(rirBright, lenInx);
    [~, inx1] = max(pks1);
%     inx1 = find(abs(max(pks1)- pks1) < 1e-6); 
    [pks2, locs2] = findpeaks(-rirBright, lenInx);
    [~, inx2] = max(pks2);

%     a = min(locs1(inx1 - 4), locs2(inx2 - 4));
%     b = max(locs1(inx1 - 4), locs2(inx2 - 4));
    space = 15;
    a = min(locs1(inx1 - space), locs2(inx2 - space));
    b = max(locs1(inx1 - space), locs2(inx2 - space));
    [~, d] = min(abs(rirBright(a:b)));
    leftInx = a + (d-1);

    a = min(locs1(inx1 + space), locs2(inx2 + space));
    b = max(locs1(inx1 + space), locs2(inx2 + space));
    [~, d] = min(abs(rirBright(a:b)));
    rightInx = a + (d-1);

%     figure()  
%     plot(leftInx, rirBright(leftInx), 'yo');hold on
%     plot(rightInx, rirBright(rightInx), 'yo');hold on
end