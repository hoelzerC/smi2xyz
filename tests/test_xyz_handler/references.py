import torch

refs = {
    "h2o.xyz": {
        "data": torch.tensor(
            [
                [8.0000, 0.0000, 0.0000, -0.7429],
                [1.0000, -1.4347, 0.0000, 0.3714],
                [1.0000, 1.4347, 0.0000, 0.3714],
            ]
        ),
        "chrg": torch.tensor([0]),
    },
    "nh3.xyz": {
        "data": torch.tensor(
            [
                [7.0000, 0.0000, 0.0000, -0.5504],
                [1.0000, 0.8849, -1.5326, 0.1835],
                [1.0000, 0.8849, 1.5326, 0.1835],
                [1.0000, -1.7697, 0.0000, 0.1835],
            ]
        ),
        "chrg": torch.tensor([0]),
    },
    "vancoh2.xyz": {
        "data": torch.tensor(
            [
                [1.00000, -4.12539, -3.21158, -3.19683],
                [6.00000, -6.02010, -2.74628, -2.59280],
                [6.00000, -7.53961, -1.20824, -4.05346],
                [6.00000, -9.93840, -0.51685, -3.22140],
                [6.00000, -10.89372, -1.63731, -1.04863],
                [6.00000, -9.38212, -3.17164, 0.43565],
                [6.00000, -6.89777, -3.61056, -0.26810],
                [8.00000, -5.30603, -4.86702, 1.38286],
                [6.00000, -3.14746, -3.62615, 2.04213],
                [6.00000, -2.82250, -1.05043, 1.59990],
                [6.00000, -0.56662, 0.11333, 2.17426],
                [6.00000, 1.31538, -1.23829, 3.38823],
                [6.00000, 0.99796, -3.78593, 3.93986],
                [6.00000, -1.22557, -5.02429, 3.20408],
                [8.00000, -1.43457, -7.58315, 3.62080],
                [6.00000, -1.56088, -9.06017, 1.35190],
                [1.00000, -2.49526, -7.95566, -0.13678],
                [6.00000, 1.14654, -9.66945, 0.46030],
                [1.00000, 2.30779, -7.97119, 0.74224],
                [8.00000, 1.00144, -10.21047, -2.15749],
                [6.00000, 3.07306, -9.30997, -3.62232],
                [6.00000, 2.42242, -9.57285, -6.42192],
                [1.00000, 1.66186, -11.45089, -6.76976],
                [1.00000, 4.17842, -9.35309, -7.47757],
                [6.00000, 0.54242, -7.57761, -7.38911],
                [7.00000, 0.58042, -7.71034, -10.13598],
                [1.00000, -0.73017, -6.51499, -10.86469],
                [1.00000, 2.29542, -7.13964, -10.77634],
                [6.00000, -2.15573, -8.18120, -6.54410],
                [1.00000, -2.24682, -8.57209, -4.53504],
                [1.00000, -2.78424, -9.85275, -7.55710],
                [1.00000, -3.42812, -6.62936, -6.99930],
                [6.00000, 1.38614, -4.91858, -6.44151],
                [8.00000, 3.46988, -4.15334, -7.89981],
                [1.00000, 3.81644, -2.39097, -7.55342],
                [1.00000, -0.19361, -3.58226, -6.72685],
                [6.00000, 1.97486, -4.96345, -3.59013],
                [6.00000, 3.12280, -2.50234, -2.63956],
                [1.00000, 4.99168, -2.28011, -3.45824],
                [1.00000, 3.32805, -2.59094, -0.60001],
                [1.00000, 1.92953, -0.90463, -3.12318],
                [1.00000, 0.18609, -5.33248, -2.58841],
                [8.00000, 3.76426, -6.86654, -2.96233],
                [1.00000, 4.75708, -10.45247, -3.17599],
                [6.00000, 2.32550, -11.81444, 2.02197],
                [1.00000, 3.99249, -12.53185, 1.01093],
                [8.00000, 3.19339, -10.93068, 4.36532],
                [1.00000, 2.15906, -9.52869, 4.92193],
                [6.00000, 0.45683, -14.04764, 2.23997],
                [8.00000, 0.30312, -15.31605, -0.08759],
                [1.00000, -0.19042, -14.10137, -1.36362],
                [1.00000, 1.21199, -15.41482, 3.59648],
                [6.00000, -2.17657, -13.20386, 3.12281],
                [6.00000, -2.38093, -12.66151, 5.97176],
                [8.00000, -4.81938, -11.83065, 6.59037],
                [1.00000, -5.23254, -10.48047, 5.43062],
                [1.00000, -2.05724, -14.40723, 7.02290],
                [1.00000, -0.97340, -11.25762, 6.55634],
                [8.00000, -3.10606, -11.12044, 1.66973],
                [1.00000, -3.52296, -14.70895, 2.65603],
                [8.00000, 2.82422, -5.20194, 5.07592],
                [6.00000, 4.99972, -3.84441, 5.65166],
                [6.00000, 5.07947, -2.38420, 7.83055],
                [6.00000, 7.06203, -0.72227, 8.18401],
                [6.00000, 8.96842, -0.49641, 6.38395],
                [6.00000, 8.96387, -2.11830, 4.31539],
                [6.00000, 6.98344, -3.78471, 3.93887],
                [1.00000, 6.92311, -4.97678, 2.28505],
                [1.00000, 10.45968, -2.02788, 2.93398],
                [6.00000, 10.91871, 1.58819, 6.68370],
                [8.00000, 13.41327, 0.78511, 6.41693],
                [1.00000, 13.66235, 0.28105, 4.66949],
                [1.00000, 10.78475, 2.31438, 8.62151],
                [6.00000, 10.31488, 3.83228, 4.86368],
                [7.00000, 7.62552, 4.38751, 4.80083],
                [6.00000, 6.74860, 6.71827, 4.29237],
                [8.00000, 8.09796, 8.57745, 4.20353],
                [6.00000, 3.96262, 6.91914, 3.50684],
                [6.00000, 4.10051, 6.80104, 0.64805],
                [6.00000, 5.54786, 4.98757, -0.58943],
                [6.00000, 6.07543, 5.15008, -3.15069],
                [6.00000, 5.01989, 7.12690, -4.56531],
                [6.00000, 3.47606, 8.88098, -3.35112],
                [6.00000, 3.05564, 8.74171, -0.78113],
                [1.00000, 1.92765, 10.17758, 0.13306],
                [1.00000, 2.70731, 10.40427, -4.46842],
                [8.00000, 5.46102, 7.44376, -7.04659],
                [1.00000, 6.80650, 6.34197, -7.61671],
                [6.00000, 7.82745, 3.28601, -4.29701],
                [6.00000, 10.38691, 3.28884, -3.63359],
                [6.00000, 12.01110, 1.49468, -4.61925],
                [6.00000, 11.12613, -0.32744, -6.29690],
                [6.00000, 8.60979, -0.30196, -7.02649],
                [6.00000, 6.95868, 1.48214, -6.03315],
                [8.00000, 4.52503, 1.37952, -6.82894],
                [1.00000, 3.40572, 2.41872, -5.80427],
                [1.00000, 7.93283, -1.69167, -8.35133],
                [8.00000, 12.64973, -2.12016, -7.28184],
                [1.00000, 14.33238, -1.97055, -6.59384],
                [1.00000, 13.97891, 1.50259, -4.07693],
                [6.00000, 11.49248, 5.40682, -2.03861],
                [7.00000, 10.80738, 5.36858, 0.61407],
                [6.00000, 11.29104, 3.40005, 2.13660],
                [8.00000, 12.35970, 1.48811, 1.43956],
                [1.00000, 10.16391, 7.01480, 1.35015],
                [6.00000, 10.76976, 7.91727, -3.24710],
                [8.00000, 10.68703, 8.24274, -5.49167],
                [8.00000, 10.30312, 9.73382, -1.56954],
                [1.00000, 9.89746, 11.31451, -2.40577],
                [1.00000, 13.56388, 5.26583, -2.19508],
                [1.00000, 6.32212, 3.39874, 0.42914],
                [7.00000, 2.29339, 5.16051, 4.83057],
                [6.00000, 0.00563, 4.25030, 4.02320],
                [8.00000, -1.82230, 4.26248, 5.39395],
                [6.00000, -0.03124, 2.84904, 1.49457],
                [7.00000, -1.85605, 3.90958, -0.23666],
                [6.00000, -1.38398, 4.01605, -2.71867],
                [8.00000, 0.57514, 3.25418, -3.67277],
                [6.00000, -3.44948, 5.21240, -4.34780],
                [7.00000, -5.92544, 4.19013, -3.77790],
                [6.00000, -8.05264, 5.08119, -4.89177],
                [8.00000, -8.01391, 6.54957, -6.64870],
                [6.00000, -10.49634, 4.19837, -3.64992],
                [7.00000, -10.15225, 4.35782, -0.92462],
                [6.00000, -11.92939, 3.92147, 0.85150],
                [8.00000, -14.00955, 3.04456, 0.39508],
                [6.00000, -11.03217, 4.54541, 3.53841],
                [7.00000, -10.33544, 7.17418, 3.87727],
                [6.00000, -12.30525, 8.95339, 3.22029],
                [1.00000, -13.98428, 8.51952, 4.33128],
                [1.00000, -11.67144, 10.85013, 3.70946],
                [1.00000, -12.84346, 8.92980, 1.21130],
                [1.00000, -8.71004, 7.57394, 2.93798],
                [6.00000, -8.82115, 2.78296, 4.20955],
                [1.00000, -7.05383, 3.63671, 3.57132],
                [1.00000, -9.08382, 0.99807, 3.20292],
                [6.00000, -8.57327, 2.21371, 7.03766],
                [6.00000, -7.81075, 4.57317, 8.50474],
                [1.00000, -9.12710, 6.10424, 8.14350],
                [1.00000, -7.77320, 4.17974, 10.52132],
                [1.00000, -5.93339, 5.16473, 7.91778],
                [6.00000, -6.60525, 0.13713, 7.42840],
                [1.00000, -4.76881, 0.78195, 6.76792],
                [1.00000, -6.44965, -0.33153, 9.42234],
                [1.00000, -7.10930, -1.57187, 6.40263],
                [1.00000, -10.40339, 1.53926, 7.73189],
                [1.00000, -12.63983, 4.16203, 4.78912],
                [1.00000, -8.60317, 5.32279, -0.34644],
                [1.00000, -11.99164, 5.50632, -4.26443],
                [6.00000, -11.38268, 1.52666, -4.60879],
                [8.00000, -13.99285, 1.31255, -4.45046],
                [1.00000, -14.50601, 1.78232, -2.73132],
                [1.00000, -10.93047, 1.43914, -6.63181],
                [1.00000, -6.07015, 2.82649, -2.45602],
                [6.00000, -3.43722, 8.11313, -4.03348],
                [6.00000, -4.79573, 9.02568, -1.68528],
                [8.00000, -5.56732, 7.64016, -0.00276],
                [7.00000, -5.10642, 11.53534, -1.57001],
                [1.00000, -6.05978, 12.28633, -0.10151],
                [1.00000, -4.59948, 12.67794, -2.99941],
                [1.00000, -1.48584, 8.78122, -4.00796],
                [1.00000, -4.39616, 8.93017, -5.66845],
                [1.00000, -2.96530, 4.80399, -6.32475],
                [1.00000, -3.41031, 4.78396, 0.48223],
                [1.00000, 1.80971, 2.92644, 0.56468],
                [1.00000, 2.16486, 5.67271, 6.67990],
                [1.00000, 3.37534, 8.85500, 3.99001],
                [1.00000, 6.41875, 2.91053, 4.88645],
                [1.00000, 11.27581, 5.52912, 5.56900],
                [1.00000, 7.09907, 0.45579, 9.85030],
                [1.00000, 3.54947, -2.50032, 9.17336],
                [1.00000, 3.06243, -0.30103, 3.85607],
                [1.00000, -4.34337, 0.02084, 0.77714],
                [1.00000, -10.08639, -3.97269, 2.17249],
                [1.00000, -12.80871, -1.26068, -0.46347],
                [1.00000, -6.84435, -0.50808, -5.84054],
            ]
        ),
        "chrg": torch.tensor([0]),
    },
}
