function mLibrary = mLibrary()
    %edit this function to modify material library, by adding branches
    %to go from MPa to psi, multiply by 145.038, 0.062428 to go to lbm/ft^3
    % material: (name, yield [psf], young [psf], density [lbm/ft^3])
    
    % Aluminum
    mLibrary.Al.a = material('2014', 96.5*145.038*144, 73*145.038*144, 2800*0.062428, 1);
    mLibrary.Al.b = material('2018', 317.1*145.038*144, 74*145.038*144, 2800*0.062428, 2);
    mLibrary.Al.c = material('2024', 75.8*145.038*144, 73*145.038*144, 2800*0.062428, 3);
    mLibrary.Al.d = material('2618-T61(SS)', 371.9*145.038*144, 74.5*145.038*144, 2760*0.062428, 4);
    mLibrary.Al.e = material('6061-T6(SS)', 275*145.038*144, 69*145.038*144, 2700*0.062428, 5);
    mLibrary.Al.f = material('6063-T83', 240*145.038*144, 69*145.038*144, 2700*0.062428, 6);
    mLibrary.Al.g = material('7050-T7651', 490*145.038*144, 72*145.038*144, 2830*0.062428, 7);
    mLibrary.Al.h = material('7075-T6(SN)', 505*145.038*144, 72*145.038*144, 2810*0.062428, 8);

    % Steel
    mLibrary.St.a = material('Alloy Steel', 620*145.038*144, 210*145.038*144, 7700*0.062428, 9);
    mLibrary.St.b = material('AISI 4340 normalized', 710*145.038*144, 205*145.038*144, 7850*0.062428, 10);

    % Titanium
    mLibrary.Ti.a = material('Ti-8Al-1Mo-1V annealed', 930.79*145.038*144, 120*145.038*144, 4369.9*0.062428, 11);

    % Fiberglass
    mLibrary.F.a = material('Amalga Fiberglass', 27000*144, 106*1.3*144, 0.072*12^3, 12);

    % Carbon fiber
    % the yield and youngs modulus are in ksi but converted to psi
    mLibrary.C.a = material('Carbon Fiber 45% resin', 515.54*1000*144, 18898*1000*144, 1600*0.062428, 13);
end
