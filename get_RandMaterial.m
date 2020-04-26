function [material] = get_RandMaterial(n)
% n is number of materials we want to sample, probably 11 or 13 (see mLibrary.m)

mat_lib = mLibrary();
material_no = round(0.5 + rand*n);

switch material_no
    case 1
        material = mat_lib.Al.a;
    case 2
        material = mat_lib.Al.b;
    case 3
        material = mat_lib.Al.c;
    case 4
        material = mat_lib.Al.d;
    case 5
        material = mat_lib.Al.e;
    case 6
        material = mat_lib.Al.f;
    case 7
        material = mat_lib.Al.g;
    case 8
        material = mat_lib.Al.h;
    case 9
        material = mat_lib.St.a;
    case 10
        material = mat_lib.St.b;
    case 11
        material = mat_lib.Ti.a;
    case 12
        material = mat_lib.F.a;
    case 13
        material = mat_lib.C.a;
end