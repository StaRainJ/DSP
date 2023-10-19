function WriteMatches(X, Y, name)
N = size(X,1);
name = [name,'.txt'];
fid = fopen(name,'a');
fprintf(fid,'%d\n',N);
for i = 1:N
    fprintf(fid,'%f %f %f %f\n',X(i,1),X(i,2),Y(i,1),Y(i,2));
end
fclose(fid);