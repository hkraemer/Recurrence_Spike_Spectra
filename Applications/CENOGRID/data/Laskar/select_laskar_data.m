clear, clc

ecc = load("La2010d_ecc3.dat");

eccentricity = ecc(1:13421,2);

save("ecc_Laskar.txt", "eccentricity", "-ascii")

