function [bif,root]=root_algorithm(equation,range)

solution=equation(range);

if solution(1)>0
    root=find(solution<=0,1); 
    bif=range(root);
else
    root=find(solution>=0,1); 
    bif=range(root);
end

