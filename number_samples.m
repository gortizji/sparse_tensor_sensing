function [s] = number_samples(select)
    s = 1;
    for i = 1:length(select)
        s = s * length(select{i});
    end
end