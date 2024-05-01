clear; clc; close all

t = [1 2 6 7 10 13 20 25]; % times (s)
v = t.^2; % corresponding velocities (ft/s)

tinstant = 8; % time for calculating instantaneous velocity
t1 = 3.0; t2 = 14; % subset of time
% alter the above lines of code

% Code to ensure tinstant, t1, and t2 are allowed
if( tinstant < t(1) | tinstant > t(end) )
    disp('WARNING: tinstant outside allowed range. Program terminating.')
    return
elseif(t1 < t(1) | t2 > t(end))
    disp('WARNING: t1 or t2 outside allowed range. Program terminating.')
    return
elseif(t1 > t2)
    disp('WARNING: t1 should be less than t2. Program terminating.')
    return
end

[totaldist,vinstant,subdist] = quadspline(t,v,tinstant,t1,t2);

disp(['Total distance traveled is ', num2str(totaldist), 'ft'])
disp(['Distance traveled from ', num2str(t1), ' to ', num2str(t2), ...
    ' seconds is ', num2str(subdist), 'ft' ])

disp(['Velocity at ', num2str(tinstant), 's is ', num2str(vinstant), 'ft/s'])