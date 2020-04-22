function mLibrary = mLibrary()
    %edit this function to modify material library, by adding branches
    %to go from MPa to psi, multiply by 145.038, 0.062428 to go to lbm/ft^3

    % Aluminum
    mLibrary.Al.one = material('2014',96.5*145.038,73*145.038,2800*0.062428);
    mLibrary.Al.two = material('2018',317.1*145.038,74*145.038,2800*0.062428);
    mLibrary.Al.three = material('2024',75.8*145.038,73*145.038,2800*0.062428);
    mLibrary.Al.four = material('2618-T61(SS)',371.9*145.038,74.4999*145.038,2760*0.062428);
    mLibrary.Al.five = material('6061-T6(SS)',275*145.038,69*145.038,2700*0.062428);
    mLibrary.Al.six = material('6063-T83',240*145.038,69*145.038,2700*0.062428);
    mLibrary.Al.seven = material('7050-T7651',490*145.038,72*145.038,2830*0.062428);
    mLibrary.Al.eight = material('7075-T6(SN) (usual choice)',505*145.038,72*145.038,2810*0.062428);

    % Steel
    mLibrary.St.one = material('7075-T6(SN) (usual choice)',505*145.038,72*145.038,2810*0.062428);
end
