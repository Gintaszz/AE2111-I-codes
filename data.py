import numpy as np
import scipy as sp
from matplotlib import pyplot as plt

SPAN = 28.15
SEMISPAN = SPAN/2
C_R = 3.13
TAPER = 0.335
DIHEDRAL = 1.3 * np.pi / 180


def c(x):
    return C_R - C_R*(1-TAPER)/(SEMISPAN) * x

n = 2.64

E = 68.9*10**9 #Pa
G = 26*10**9 #Pa

airfoil_data_top = [1.00000, 0.00000,
0.95032, 0.00976,
0.90066, 0.01982,
0.85092, 0.03020,
0.80109, 0.04062,
0.75115, 0.05084,
0.70111, 0.06055,
0.65096, 0.06954,
0.60072, 0.07762,
0.55040, 0.08456,
0.50000, 0.09016,
0.44954, 0.09414,
0.39904, 0.09614,
0.34853, 0.09541,
0.29803, 0.09260,
0.24756, 0.08771,
0.19714, 0.08066,
0.14681, 0.07122,
0.09662, 0.05864,
0.07162, 0.05075,
0.04673, 0.04121,
0.02207, 0.02883,
0.00996, 0.02038,
0.00526, 0.01579,
0.00299, 0.01291,
0.00000, 0.00000]

airfoil_data_top.reverse()
for index in np.arange(0, len(airfoil_data_top)-1, step=2):
    airfoil_data_top[index + 1], airfoil_data_top[index] = airfoil_data_top[index], airfoil_data_top[index + 1]

airfoil_data_bottom = [0.00000, 0.00000,
0.00701, -0.01091,
0.00974, -0.01299,
0.01504, -0.01610,
0.02793, -0.02139,
0.05327, -0.02857,
0.07838, -0.03379,
0.10338, -0.03796,
0.15319, -0.04430,
0.20286, -0.04882,
0.25244, -0.05191,
0.30197, -0.05372,
0.35147, -0.05421,
0.40096, -0.05330,
0.45046, -0.05034,
0.50000, -0.04604,
0.54960, -0.04076,
0.59928, -0.03478,
0.64904, -0.02834,
0.69889, -0.02167,
0.74885, -0.01504,
0.79891, -0.00878,
0.84908, -0.00328,
0.89934, 0.00086,
0.94968, 0.00288,
1.00000, 0.00000]

airfoilfunc_top = sp.interpolate.InterpolatedUnivariateSpline(
                [airfoil_data_top[i] for i in np.arange(0, len(airfoil_data_top), 2)], 
                [airfoil_data_top[i+1] for i in np.arange(0, len(airfoil_data_top), 2)]
                )

airfoilfunc_bottom = sp.interpolate.InterpolatedUnivariateSpline(
                [airfoil_data_bottom[i] for i in np.arange(0, len(airfoil_data_bottom), 2)], 
                [airfoil_data_bottom[i+1] for i in np.arange(0, len(airfoil_data_bottom), 2)]
                )

designparameters = {'area stringer': [2*10**-4, 2*10**-4, 2*10**-4], 
                    't spar': [0.01, 0.02, 0.02], 
                    't web': [0.001, 0.01, 0.0005], 
                    'front spar x': [0.25, 0.25, 0.25], 
                    'back spar x': [0.67, 0.67, 0.67], 
                    'list stringers': [[24, 18, 14, 10, 6, 0],
                                       [24, 14, 10, 6, 0],
                                       [24, 18, 14, 10, 6, 0]],
              }

designproperties = {'span list': [],
                    'spar distance x': [],
                    'front spar to root chord': [],
                    'back spar to root chord': [],
                    }

#Create design parameters dictionary
for designindex in range(3):
    designproperties['span list'].append(\
        [i*SEMISPAN/(len(designparameters['list stringers'][designindex])-1)\
        for i in range(len(designparameters['list stringers'][designindex]))])
    
    designproperties['spar distance x'].append(designparameters['back spar x'][designindex] - designparameters['front spar x'][designindex])
    
    designproperties['front spar to root chord'].append((airfoilfunc_top(designparameters['front spar x'][designindex]) - 
                                                         airfoilfunc_bottom(designparameters['front spar x'][designindex])) * C_R)
    
    designproperties['back spar to root chord'].append((airfoilfunc_top(designparameters['back spar x'][designindex]) - 
                                                         airfoilfunc_bottom(designparameters['back spar x'][designindex])) * C_R)

#Merge dictionaries into 1 named design
design = {**designparameters, **designproperties}

def topweb(x, c, designnum): #Equation of top web from geometry chord as datum
    frontx, fronty = design['front spar x'][designnum], airfoilfunc_top(design['front spar x'][designnum])
    backx, backy = design['back spar x'][designnum], airfoilfunc_top(design['back spar x'][designnum])
    slope = (fronty-backy)/(frontx-backx)
    return slope*x + fronty*c


def bottomweb(x, c, designnum): #Equation of top web from geometry chord as datum
    frontx, fronty = design['front spar x'][designnum], airfoilfunc_bottom(design['front spar x'][designnum])
    backx, backy = design['back spar x'][designnum], airfoilfunc_bottom(design['back spar x'][designnum])
    slope = (fronty-backy)/(frontx-backx)
    return slope*x + fronty*c


def centroid(c, designnum):
    t_spar = design['t spar'][designnum]
    t_web = design['t web'][designnum]
    spardist_over_chord = design['back spar x'][designnum] - design['front spar x'][designnum]

    centroidfrontsparcr = (airfoilfunc_top(design['front spar x'][designnum]) + 
                           airfoilfunc_bottom(design['front spar x'][designnum])) /2
    
    centroidbacksparcr = (airfoilfunc_top(design['back spar x'][designnum]) + 
                          airfoilfunc_bottom(design['back spar x'][designnum])) /2

    total_area = design['front spar to root chord'][designnum]/C_R * t_spar + \
                design['back spar to root chord'][designnum]/C_R * t_spar + \
                spardist_over_chord * t_web + \
                spardist_over_chord * t_web

    frontsparcontr = design['front spar to root chord'][designnum]/C_R * t_spar * centroidfrontsparcr
    backsparcontr = design['back spar to root chord'][designnum]/C_R * t_spar * centroidbacksparcr
    topwebcentcontr = spardist_over_chord * t_web * (topweb(0, 1, designnum) + topweb(spardist_over_chord, 1, designnum))/2
    bottomwebcentcontr = spardist_over_chord * t_web * (bottomweb(0, 1, designnum) + bottomweb(spardist_over_chord, 1, designnum))/2    
    cen = (frontsparcontr + backsparcontr + topwebcentcontr + bottomwebcentcontr)/(total_area)
    return cen * c


if __name__ == '__main__':
    #-----------------------------Plot winwbox unit chord---------------------------------------
    res = 100

    fig, axs = plt.subplots(3)

    for designindex in range(len(design['front spar x'])):
        axs[designindex].set_title(f'Design {designindex+1}')
        axs[designindex].plot(np.arange(0, 1, 0.001), [airfoilfunc_top(x) for x in np.arange(0, 1, 0.001)], color = 'r')
        axs[designindex].plot(np.arange(0, 1, 0.001), [airfoilfunc_bottom(x) for x in np.arange(0, 1, 0.001)], color = 'r')

        sampling = np.linspace(design['front spar x'][designindex], design['back spar x'][designindex], res)
        axs[designindex].plot(sampling, [topweb(x - design['front spar x'][designindex], 1, designindex) for x in sampling], color = 'b')
        axs[designindex].plot(sampling, [bottomweb(x - design['front spar x'][designindex], 1, designindex) for x in sampling], color = 'b')
        axs[designindex].vlines(design['front spar x'][designindex], ymin = bottomweb(0, 1, designindex), ymax = topweb(0, 1, designindex), color = 'b')
        axs[designindex].vlines(design['back spar x'][designindex], ymin = bottomweb(design['spar distance x'][designindex], 1, designindex), ymax = topweb(design['spar distance x'][designindex], 1, designindex), color = 'b')
        axs[designindex].axhline(centroid(1, designindex), color = 'g')
        axs[designindex].axis('equal')

    plt.show()

#-------------------------------------------------------------------------------------------------

momentwing = np.array([-262354.8292325158, -261740.76316955494, -260902.0549120349, -259839.1247890445, -258552.39573383873, -257042.2807733308, -255309.22891908276, -253353.5504149055, -251176.07428002107, -248774.94003372442, -246152.52629418386, -243311.16949213253, -240249.1090321707, -236966.7644608659, -233465.89750047115, -229750.61847486257, -225956.27682770224, -222194.32121234684, -218448.84569974057, -214724.57043137212, -211020.65800397593, -207337.80815979448, -203676.06367084553, -200036.49153336257, -196416.65688049115, -192827.04714323225, -189229.9065789975, -185855.77944033738, -182704.93344768533, -179547.36194801267, -176421.3525617724, -173316.96249452163, -170237.16406487848, -167181.45863977977, -164150.27189930395, -161143.772425591, -158162.19409257636, -155205.74316872182, -152274.6265808048, -149369.0434808458, -146489.18932705285, -143635.2522886911, -140807.41499490794, -138005.85460175519, -135230.74257829835, -132482.24457813762, -129760.52096276495, -127065.72647782939, -124398.01041862361, -121757.51669306749, -119144.38417594018, -116558.74659178033, -114000.73270941875, -111470.46614413729, -108968.06587470809, -106493.64593726369, -104047.3158098531, -101629.18036361343, -99239.33979937343, -96877.8897933735, -94544.9213821393, -92240.5213651292, -89964.77194079834, -87717.75080130705, -85499.53135240519, -83310.18268157986, -81149.76951326283, -79018.35228915229, -76915.98709841474, -74842.72601406198, -72798.61668671732, -70783.70272387468, -68798.02355349161, -66841.6144585647, -64914.50669074504, -63016.72725301616, -61148.29935107518, -59309.24200496001, -57499.57022130769, -55719.29514154388, -53968.423910005564, -52246.95962036415, -50554.90146706076, -48892.244720755225, -47258.98066310285, -45655.09664379044, -44080.57616611548, -42535.398745840175, -41019.539948390906, -39532.976038415116, -38075.66136006687, -36647.573146829745, -35248.66689912744, -33878.898929764895, -32538.220066289465, -31226.5794650774, -29943.920711763287, -28690.18373890878, -27465.304603873457, -26269.21508883042, -25101.84300365566, -23963.111948304584, -22852.941282505708, -21771.246141461266, -20717.937305848052, -19692.921168640813, -18696.09974149597, -17727.370512845093, -16786.626426229377, -15873.75571580828, -14988.641840707423, -14131.16332277847, -13301.193635361264, -12498.601125328225, -11723.248806258141, -10974.994348111404, -10253.689836233454, -9559.181376349434, -8891.308983373525, -8249.905965136639, -7634.798863695017, -7045.807290254192, -6482.741870320934, -5945.40744201544, -5433.599083230827, -4947.103128551965, -4485.69639388422, -4049.1455113316465, -3637.2063233639333, -3249.6233010227993, -2886.1285591283777, -2546.439820959796, -2230.257800589951, -1937.2635887788272, -1667.1160269829638, -1419.4490621024245, -1193.8691846931888, -989.9529737689714, -807.2446364214611, -645.2535384484596, -503.4516759404901, -381.2621649905378, -278.0309903429889, -192.9920199287056, -125.23196608219803, -73.67060911884349, -37.06449344790329, -14.020215856527155, -2.721573148764266, 0.0])
torquewing = np.array([-71266.38846420322, -70933.9400868938, -70604.29620132309, -70277.89831389672, -69955.18794803531, -69636.6066535209, -69322.59600231075, -69013.59636847896, -68710.02331745488, -68412.23448340235, -68120.52466833989, -67835.12598979913, -67556.20830291245, -67283.87853160847, -67018.18161562277, -16647.10016112757, -16394.554685554594, -16148.403680903262, -15908.443723471399, -15674.409920002712, -15445.995462375387, -15222.903527038256, -15004.86234212427, -14791.625410488778, -14582.971258408741, -14378.703507442673, -14178.650786004784, -7408.848625080101, -7216.811580098341, -7028.624916977161, -6844.2167144959685, -6663.539769899072, -6486.567083638063, -6313.274800170647, -6143.6350819489535, -5977.61594147642, -5815.1812875275955, -5656.290710077644, -5500.899894635381, -5348.960229809828, -5200.418977156351, -5055.219233634977, -4913.299995982622, -4774.596050528893, -4639.039264850843, -4506.562458610686, -4377.100948012598, -4250.592527066731, -4126.977449946326, -4006.1985603740804, -3888.2008050880168, -3772.932015854888, -3660.3422035029585, -3550.3838456367844, -3443.011844393839, -3338.183461226032, -3235.857551675068, -3135.99334172368, -3038.550230374741, -2943.487788281675, -2850.765748995618, -2760.344011641555, -2672.182627033287, -2586.241800934557, -2502.481894296207, -2420.863405557699, -2341.346983848358, -2263.8934997433166, -2188.464326410315, -2115.0213969693486, -2043.527296190828, -1973.945056745588, -1906.2383996825206, -1840.3715641418412, -1776.309363719215, -1714.01717317293, -1653.4609505318547, -1594.607172282231, -1537.422830852062, -1481.875284310007, -1427.9322180138295, -1375.5616684891702, -1324.7319734007567, -1275.411840411924, -1227.5703079531104, -1181.1767159178132, -1136.2007561264063, -1092.612448930876, -1050.382136506677, -1009.4804765974815, -969.8784026759012, -931.5471866580705, -894.4583792767012, -858.5838450350725, -823.8957523635864, -790.366574240606, -757.9691067261044, -726.6764551795197, -696.4619941561763, -667.2993665732208, -639.1624436304523, -612.0253526003876, -585.8624646750516, -560.6484094972662, -536.3580833472715, -512.9666481385894, -490.44954289226007, -468.7824363264921, -447.94110894280306, -427.9014179135273, -408.63931928941275, -390.13087656419805, -372.3522684759541, -355.2798421485573, -338.890047255676, -323.15951547830326, -308.06497327977877, -293.5832318491428, -279.69122532687953, -266.3660501946566, -253.5850014157312, -241.32561501025486, -229.56569348476347, -218.2829817803592, -207.45440040127176, -197.05591641625878, -187.0626098730921, -177.4487467978476, -168.18783988786365, -159.2527926083743, -150.6169564622086, -142.25559234135147, -134.14631756839654, -126.26941122467923, -118.60811732380861, -111.14799538238151, -103.8674849722048, -96.72896696992967, -89.67826912738695, -82.64494595719232, -75.54743420962605, -68.33720465253455, -61.03850213207677, -53.75502769907957, -46.639474705899715, -39.58831284121659, -31.978566863369217, -22.763105712928663, -11.608696872725336, 0.0])

moment = sp.interpolate.InterpolatedUnivariateSpline(np.linspace(0, SEMISPAN, len(momentwing)), n*momentwing)
torque = sp.interpolate.InterpolatedUnivariateSpline(np.linspace(0, SEMISPAN, len(torquewing)), n*torquewing, k=1)