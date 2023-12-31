## table of equivalent reduced coordinate sets for the 38 canonical puckers
##
##       "Berce's et al. dihedrals" "Cremer-Pople coordinates" "puckering angles"
## conformer tau1 tau2 tau3              phi theta BigO       theta1 theta2 theta3
berces = [
    ["1C4", 60, -60, 60],
    ["4C1", -60, 60, -60],
    ["1,4B", 0, 60, -60],
    ["B1,4", 0, -60, 60],
    ["2,5B", -60, 0, 60],
    ["B2,5", 60, 0, -60],
    ["3,6B", 60, -60, 0],
    ["B3,6", -60, 60, 0],
    ["1H2", 45, -15, 0],
    ["2H1", -45, 15, 0],
    ["2H3", -60, 45, -15],
    ["3H2", 60, -45, 15],
    ["3H4", 45, -60, 45],
    ["4H3", -45, 60, -45],
    ["4H5", -15, 45, -60],
    ["5H4", 15, -45, 60],
    ["5H6", 0, -15, 45],
    ["6H5", 0, 15, -45],
    ["6H1", -15, 0, -15],
    ["1H6", 15, 0, 15],
    ["1S3", -30, 60, -30],
    ["3S1", 30, -60, 30],
    ["5S1", -30, -30, 60],
    ["1S5", 30, 30, -60],
    ["6S2", 60, -30, -30],
    ["2S6", -60, 30, 30],
    ["1E", 30, 0, 0],
    ["E1", -30, 0, 0],
    ["2E", -60, 30, 0],
    ["E2", 60, -30, 0],
    ["3E", 60, -60, 30],
    ["E3", 60, 60, -30],
    ["4E", -30, 60, -60],
    ["E4", 30, -60, 60],
    ["5E", 0, -30, 60],
    ["E5", 0, 30, -60],
    ["6E", 0, 0, -30],
    ["E6", 0, 0, 30],
]

KremerPol = [
    ["1C4", 0, 180, 0.57],
    ["4C1", 0, 0, 0.57],
    ["1,4B", 240, 90, 0.76],
    ["B1,4", 60, 90, 0.76],
    ["2,5B", 120, 90, 0.76],
    ["B2,5", 300, 90, 0.76],
    ["3,6B", 0, 90, 0.76],
    ["B3,6", 180, 90, 0.76],
    ["1H2", 270, 129, 0.42],
    ["2H1", 90, 51, 0.42],
    ["2H3", 150, 51, 0.42],
    ["3H2", 330, 129, 0.42],
    ["3H4", 30, 129, 0.42],
    ["4H3", 210, 51, 0.42],
    ["4H5", 270, 51, 0.42],
    ["5H4", 90, 129, 0.42],
    ["5H6", 150, 129, 0.42],
    ["6H5", 330, 51, 0.42],
    ["6H1", 30, 51, 0.42],
    ["1H6", 210, 129, 0.42],
    ["1S3", 210, 88, 0.62],
    ["3S1", 30, 92, 0.62],
    ["5S1", 90, 92, 0.62],
    ["1S5", 270, 88, 0.62],
    ["6S2", 330, 88, 0.62],
    ["2S6", 150, 92, 0.62],
    ["1E", 240, 125, 0.45],
    ["E1", 60, 55, 0.45],
    ["2E", 120, 55, 0.45],
    ["E2", 300, 125, 0.45],
    ["3E", 360, 125, 0.45],
    ["E3", 180, 55, 0.45],
    ["4E", 240, 55, 0.45],
    ["E4", 60, 125, 0.45],
    ["5E", 120, 125, 0.45],
    ["E5", 300, 55, 0.45],
    ["6E", 360, 55, 0.45],
    ["E6", 180, 125, 0.45],
]
HillReilly = [
    ["1C4", -35.26, -35.26, -35.26],
    ["4C1", 35.26, 35.26, 35.26],
    ["1,4B", -35.26, 74.20, -35.26],
    ["B1,4", 35.26, -74.20, 35.26],
    ["2,5B", 74.20, -35.26, -35.26],
    ["B2,5", -74.20, 35.26, 35.26],
    ["3,6B", -35.26, -35.26, 74.20],
    ["B3,6", 35.26, 35.26, -74.20],
    ["1H2", -42.16, 9.07, -17.83],
    ["2H1", 42.16, -9.07, 17.83],
    ["2H3", 42.16, 17.83, -9.06],
    ["3H2", -42.16, -17.83, 9.06],
    ["3H4", -17.83, -42.16, 9.07],
    ["4H3", 17.83, 42.16, -9.07],
    ["4H5", -9.07, 42.16, 17.83],
    ["5H4", 9.07, -42.16, -17.83],
    ["5H6", 9.07, -17.83, -42.16],
    ["6H5", -9.07, 17.83, 42.16],
    ["6H1", 17.83, -9.07, 42.16],
    ["1H6", -17.83, 9.07, -42.16],
    ["1S3", 0, 50.84, -50.84],
    ["3S1", 0, -50.84, 50.84],
    ["5S1", 50.84, -50.84, 0],
    ["1S5", -50.84, 50.84, 0],
    ["6S2", -50.84, 0, 50.84],
    ["2S6", 50.84, 0, -50.84],
    ["1E", -35.26, 17.37, -35.26],
    ["E1", 35.26, -17.37, 35.26],
    ["2E", 46.86, 0, 0],
    ["E2", -46.86, 0, 0],
    ["3E", -35.26, -35.26, 17.37],
    ["E3", 35.26, 35.26, -17.37],
    ["4E", 0, 46.86, 0],
    ["E4", 0, -46.86, 0],
    ["5E", 17.37, -35.26, -35.26],
    ["E5", -17.37, 35.26, 35.26],
    ["6E", 0, 0, 46.86],
    ["E6", 0, 0, -46.86],
]
HillReillyShort = [
    HillReilly[0],
    HillReilly[1],
    HillReilly[3],
    HillReilly[4],
    HillReilly[6],
    HillReilly[7],
]

HillReillyQMOpt6 = [
    ["1C4", [-32.00, -32.00, -32.00]],
    ["4C1", [32.00, 32.00, 32.00]],
]
HillReillyQMOpt5 = [["1C4", [-32.00, 32.00]], ["4C1", [32.00, -32.00]]]
## # CH04
##      0 (-0.35259309232847613, 56.948716515764865, -56.59590592826308)
##      1 (-56.56301095755276, -0.33274504416716866, 56.96756218902601)
##      2 (56.94871651576488, -56.59590592826311, -0.35259309232847613)
## # CH06
##      0 (63.62450768150437, -30.424699015618586, -31.364036922770524)
##      1 (-31.41055465170635, 63.57871950280473, -30.49607246803778)
##      2 (-30.424699015618586, -31.364036922770524, 63.62450768150437)
## # CH07
##      0 (-56.708295466381855, 56.77435164776965, -0.22882384348595508)
##      1 (-0.09093847641567265, -56.57295297561805, 56.82103681151629)
##      2 (56.77435164776965, -0.2288238434859693, -56.708295466381855)
