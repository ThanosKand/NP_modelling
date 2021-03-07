my_time = 0;
for i = 1:1000
    if round(P(i+1, end)) == round(P(i, end))
        my_time = t1(i)
%         break;
    end
end
    
    
    