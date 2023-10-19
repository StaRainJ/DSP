function d = TransferDistanceH(pts1h, pts2h, H)
V = H*pts1h;
V = V*(diag(1./V(3,:)));
d = sqrt(sum((pts2h - V).^2)');
% d = SampsonDistanceH(H,pts1h,pts2h);
end