function c = cumulative(a,b)

for i = 1:length(b)
    c(i) = length(find(a<=b(i)))/length(a);
end