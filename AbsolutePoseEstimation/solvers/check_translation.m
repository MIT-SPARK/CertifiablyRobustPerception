function flag = check_translation(t,tBound,FOV)
flag = true;
if norm(t) > tBound
    flag = false;
end
if t(3) < 0
    flag = false;
end
halfFOV = deg2rad(FOV/2);

if (tan(halfFOV))^2 * t(3)^2 - (t(1)^2 + t(2)^2) < 0
    flag = false;
end
end