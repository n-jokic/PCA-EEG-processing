function r = rmss(s)
% racuna snagu signala po formuli 
    r = 0; 
    s = s.^2; 
    r = sum(sum(s));
    r = r/size(s,1);
    r = r/size(s,2); 
    r = sqrt(r); 

