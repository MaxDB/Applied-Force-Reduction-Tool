function rotation_matrix = abaqus_rotation(rotation)
%https://help.3ds.com/2026x/english/dsdoc/SIMA3DXMODRefMap/simamod-c-conventions.htm?contextscope=cloud&id=6c63a4d7177f40fa8a6e80bce2c7d2e7#simamod-c-conventions-t-DegreesOfFreedom-sma-topic1
%abqus docs: simulation/structures:abaqus/conventions/finite rotations 
rot_angle = sqrt(sum(rotation.^2));
if rot_angle == 0
    rotation_matrix = eye(3);
    return
end
rot_axis = rotation/rot_angle;

ct = cos(rot_angle);
st = sin(rot_angle);

ux = rot_axis(1);
uy = rot_axis(2);
uz = rot_axis(3);

rotation_matrix = [
    ux^2*(1-ct) + ct,       ux*uy*(1-ct) - uz*st,   ux*uz*(1-ct) + uy*st;
    ux*uy*(1-ct) + uz*st,   uy^2*(1-ct) + ct,       uy*uz*(1-ct) - ux*st;
    ux*uz*(1-ct) - uy*st,   uy*uz*(1-ct) + ux*st,   uz^2*(1-ct) + ct
];


% u_cross = [
%     0,      -uz,    uy;
%     uz,     0,      -ux;
%     -uy,    ux,     0];
% 
% u_outer = rot_axis'.*rot_axis;
% rotation_matrix = ct*eye(3) + st*u_cross + (1-ct)*u_outer;
end