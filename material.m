function material = material(name, yield, young, density, no)
    %material(yield,young,density)
    %Input arguments to format into material struct to be stored in library
    material = struct('name',[],'yield',[],'young',[],'density',[],'no',[]);
    material.name = name;
    material.yield = yield;
    material.young = young;
    material.density = density;
    material.no = no;
end
