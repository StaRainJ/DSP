function plot_matches(I1, I2, X, Y, VFCIndex, CorrectIndex)

NumPlot = 100;

TruePos = intersect(VFCIndex, CorrectIndex);%Ture positive
FalsePos = setdiff(VFCIndex, CorrectIndex); %False positive
FalseNeg = setdiff(CorrectIndex, VFCIndex); %False negative

NumPos = length(TruePos)+length(FalsePos)+length(FalseNeg);
if NumPos > NumPlot
    t_p = length(TruePos)/NumPos;
    n1 = round(t_p*NumPlot);
    f_p = length(FalsePos)/NumPos;
    n2 = ceil(f_p*NumPlot);
    f_n = length(FalseNeg)/NumPos;
    n3 = ceil(f_n*NumPlot);
else
    n1 = length(TruePos);
    n2 = length(FalsePos);
    n3 = length(FalseNeg);
end

per = randperm(length(TruePos));
TruePos = TruePos(per(1:n1));
per = randperm(length(FalsePos));
FalsePos = FalsePos(per(1:n2));
per = randperm(length(FalseNeg));
FalseNeg = FalseNeg(per(1:n3));

%I1(size(I2,1),size(I2,2),:)=0;


[M1,N1,P1]=size(I1);
[M2,N2,P2]=size(I2);

%%
if size(I1,3) == size(I2,3)
else
if size(I1,3)==1
    I1 = repmat(I1,[1,1,3]);
end
if size(I2,3)==1
    I2 = repmat(I2,[1,1,3]);
end
end

[wa,ha,~] = size(I1);
[wb,hb,~] = size(I2);
maxw = max(wa,wb);maxh = max(ha,hb);
I2(wb+1:maxw, :,:) = 0;
I1(wa+1:maxw, :,:) = 0;

%%
% I2=imresize(I2,[M1,N1]);
% Y(:,1)=Y(:,1)*M1/M2;
% Y(:,2)=Y(:,2)*N1/N2;

interval = 10;
WhiteInterval = 255*ones(size(I1,1), interval, size(I1,3));
% figure;imagesc(cat(2, I1, WhiteInterval, I2)) ;
imshow(cat(2, I1, WhiteInterval, I2),[],'border','tight')


hold on ;
line([X(FalsePos,1)'; Y(FalsePos,1)'+size(I1,2)+interval], [X(FalsePos,2)' ;  Y(FalsePos,2)'],'linewidth', 1.5, 'color', [0.9,0.05,0.05]) ;
line([X(FalseNeg,1)'; Y(FalseNeg,1)'+size(I1,2)+interval], [X(FalseNeg,2)' ;  Y(FalseNeg,2)'],'linewidth', 1.5, 'color', [0.05,0.9,0.05]) ;
line([X(TruePos,1)'; Y(TruePos,1)'+size(I1,2)+interval], [X(TruePos,2)' ;  Y(TruePos,2)'],'linewidth', 1.5, 'color', [0.05,0.05,0.9]) ;
axis equal ;axis off  ; 
drawnow;

cols = size(I1,2);

  ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];

set(gca,'position',[0 0 1 1])
axis([0 N1+cols+interval 0 M1])
axis equal ;axis off  ; 
drawnow;








