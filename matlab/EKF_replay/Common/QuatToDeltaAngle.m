% Convert from a quaternion to a delta angle in radians

function deltaAngle = QuatToDeltaAngle(quat)

if (quat(1) >= 0.0)
    delta = 2 * acos(quat(1));
    x = quat(2) / sin(delta/2);
    y = quat(3) / sin(delta/2);
    z = quat(4) / sin(delta/2);
else
    delta = 2 * acos(-quat(1));
    x = -quat(2) / sin(delta/2);
    y = -quat(3) / sin(delta/2);
    z = -quat(4) / sin(delta/2);
end

xyz_norm = sqrt(x*x + y*y + z*z);
deltaAngle = delta / xyz_norm * [x;y;z];

