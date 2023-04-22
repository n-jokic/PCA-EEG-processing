function [ss, yy] = matching_components(s,y)
% function returns indexes of y correspondig to s
    n = size(y,1);
    m = size(s,1);
%
    for i = 1:n
        y(i,:) = ((y(i,:)-mean(y(i,:)))./(std(y(i,:))));
    end
%
    temp = zeros(n,m); 
    for i = 1:n
        for j = 1:m
            temp(i,j) = abs(corr(s(j,:)',y(i,:)'));
        end
    end    
    %    
    [pom , index] = max(temp);
    %
    yy = zeros(size(s,1),size(y,2)); 
    %
    for i = 1:size(s,1); 
        yy(i,:) = y(index(i),:);
    end
    %
    ss = s; 