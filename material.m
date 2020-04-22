function material = material(name,yield,young,density)
    %material(yield,young,density)
    %Input arguments to format into material struct to be stored in library
    material = struct('name',[],'yield',[],'young',[],'density',[]);
    material.name = name;
    material.yield = yield;
    material.young = young;
    material.density = density;
end
