function [Tmid] = FUN_get_T_from_h(h)

% this function maps T back from h

max_iter = 2000;
Thigh = 6000;
Tlow = 200;
res_h = 100;
iter1 = 0;

while res_h > 10^-4
    iter1=iter1+1;
    Tmid = (Thigh+Tlow)/2;
    hmid = FUN_get_prop_NASA9(Tmid,'h');
    if hmid>h
        Thigh = Tmid;
    else
        Tlow = Tmid;
    end
    res_h = sqrt(((h-hmid)/h)^2);
    if iter1 > max_iter
        break;
    end
end

end