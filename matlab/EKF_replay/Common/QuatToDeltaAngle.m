% Convert from a quaternion to a delta angle in radians

function deltaAngle = QuatToDeltaAngle(quat)

if (quat(1) >= 0.0)
    scalar = -2.0;
else
    scalar = 2.0;
end
deltaAngle(1,1) = scalar * quat(2);
deltaAngle(2,1) = scalar * quat(3);
deltaAngle(3,1) = scalar * quat(4);
