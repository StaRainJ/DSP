function [H, degeneracy] = DegeneracyCheck(pts1h, pts2h, indices, th)

degenIndices = [0, 1, 2, 3; 3, 4, 5, 6; 0, 1, 5, 6; 0, 2, 4, 5; 1, 2, 4, 6; 0, 3, 4, 6; 1, 3, 4, 5; 2, 3, 5, 6];
testIndices = [4, 5, 6; 0, 1, 2; 2, 3, 4; 1, 3, 6; 0, 3, 5; 1, 2, 5; 0, 2, 6; 0, 1, 4];

degeneracy = false;
for i = 1:8
    H =  norm4Point(pts1h(:,indices(degenIndices(i,:)+1)),pts2h(:,indices(degenIndices(i,:)+1)));
    d = SampsonDistanceH(pts1h(:,indices(testIndices(i,:)+1)),pts2h(:,indices(testIndices(i,:)+1)),H);
    if ~all(d>2.0*th)
        degeneracy = true;
        break;
    end
end

    