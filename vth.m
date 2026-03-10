function v=vth(T)
    mass = 9.10938215e-31; % Same as Gkeyll
    e=1.60217487e-19;
    v = sqrt(2*(T*e)/mass);
end 