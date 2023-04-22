function r = rmse(s,y)
    p = size(s,1);
    r1 = zeros(p,1); 
    r2 = zeros(p,1); 
    y_inv = -y; 
for i = 1:p 
    for k = 1:size(s,2)
        r1(i) = r1(i) + ( s(i,k)- y(i,k) )^2; 
        r2(i) = r2(i) + ( s(i,k)- y_inv(i,k) )^2;
    end
        r(i) = min(r1(i),r2(i));
        % in case that signal is "inverted"
        r(i) = r(i)/size(s,2); 
        r(i) = sqrt(r(i)); 
end
% DELETE INVERTED SIGNAL AND USE THIS Fnct.