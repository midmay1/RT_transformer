# RT_transformer
transformer through atom environment can be good for RT prediction?


1. Pre-processing
 - InChI-to-smiles converting (80038 --> 79957, 81 mismatches)
 - CID : [312, 325, 367, 418, 703, 808, 966, 1017, 1040, 1095, 1230, 1277, 1318, 1540, 1817, 2072, 2095, 2096, 2127, 2193, 2207, 2333, 2423, 2465, 2503, 2507, 2511, 2596, 2651, 2700, 2784, 2888, 2977, 3131, 3157, 3228, 3287, 3314, 3628, 3673, 3705, 4084, 4299, 4628, 4871, 4996, 5299, 5387, 5688, 6051, 6086, 6248, 6647, 7062, 7063, 7075, 7150, 7167, 7209, 7374, 7376, 7448, 7451, 7469, 7517, 14869, 16047, 16151, 16165, 16221, 16350, 16358, 16414, 16440, 18590, 22783, 26973, 31505, 32434, 33140, 63505]


 - smiles-to-SMARTS converting (add Hydrogen, assign ring or not, the number of closest neighbor(need to make sentence more clearly), and aromaticity??)
