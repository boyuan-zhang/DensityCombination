function x = make_num2str(x)
% apply num2str if the element is not str)

[d1,d2] = size(x);
for i=1:1:d1
    for j=1:1:d2
        
        if ~isstr(x{i,j})
            x{i,j} = num2str(x{i,j});
        end
        
    end
end
