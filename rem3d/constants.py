# python 3 compatibility
from __future__ import absolute_import, division, print_function
from builtins import *

import scipy.constants
import ConfigParser
import pkgutil
import StringIO

"""
download files
"""
downloadpage = 'https://maurya.umd.edu/files'
localfilefolder = 'files'

"""
Mapping constants
"""
geoco=0.993277    #correction for geographic-geocentric conversion


"""
Color scales
"""
colorscale = {
'rem3d': {'name': 'rem3d','description': \
	'REM3D color scale created by combining two perceptually uniform \
	 color scales linear_kry_5_95_c72 and linear_bgy_10_95_c74 from \
	 https://bokeh.github.io/colorcet',\
'RGB': [ (0, 0.047803, 0.4883),(0, 0.049756, 0.49556),(0, 0.051605, 0.50281),
(0, 0.053545, 0.51004),(0, 0.055585, 0.51721),(0, 0.057566, 0.52435),
(0, 0.05978, 0.53144),(0, 0.061812, 0.53849),(0, 0.064016, 0.5455),
(0, 0.066232, 0.55245),(0, 0.068551, 0.55934),(0, 0.070824, 0.5662),
(0, 0.073229, 0.57299),(0, 0.07557, 0.57971),(0, 0.078003, 0.58638),
(0, 0.080554, 0.59299),(0, 0.083114, 0.59951),(0, 0.085697, 0.60596),
(0, 0.08842, 0.61236),(0, 0.091168, 0.61866),(0, 0.093925, 0.62487),
(0, 0.096707, 0.63101),(0, 0.09963, 0.63705),(0, 0.1026, 0.64298),
(0, 0.10561, 0.64881),(0, 0.10866, 0.65454),(0, 0.11183, 0.66016),
(0, 0.11497, 0.66567),(0, 0.11829, 0.67103),(0, 0.12156, 0.67626),
(0, 0.12498, 0.68134),(0, 0.12846, 0.68629),(0, 0.13201, 0.69107),
(0, 0.13559, 0.6957),(0, 0.13927, 0.70014),(0, 0.14307, 0.70439),
(0, 0.1469, 0.70845),(0, 0.15085, 0.71227),(0, 0.15487, 0.71588),
(0, 0.159, 0.71923),(0, 0.16323, 0.7223),(0, 0.16754, 0.7251),
(0, 0.17195, 0.72757),(0, 0.17647, 0.72972),(0, 0.18113, 0.73149),
(0, 0.18594, 0.73289),(0, 0.19074, 0.73398),(0, 0.19556, 0.73486),
(0, 0.20033, 0.73556),(0, 0.20512, 0.73608),(0, 0.20987, 0.73643),
(0, 0.21461, 0.73659),(0, 0.21934, 0.73657),(0, 0.22402, 0.73637),
(0, 0.22875, 0.73599),(0, 0.2334, 0.73544),(0, 0.23809, 0.73469),
(0, 0.24275, 0.73376),(0, 0.24743, 0.73266),(0, 0.25208, 0.73137),
(0, 0.25673, 0.72991),(0, 0.26137, 0.72825),(0, 0.26603, 0.72642),
(0, 0.27068, 0.72441),(0, 0.27531, 0.72221),(0, 0.27995, 0.71983),
(0, 0.28458, 0.71727),(0, 0.28924, 0.71452),(0, 0.29387, 0.71161),
(0, 0.29852, 0.70851),(0, 0.30317, 0.70521),(0, 0.30782, 0.70174),
(0, 0.31248, 0.69809),(0, 0.31716, 0.69426),(0, 0.32182, 0.69025),
(0, 0.32649, 0.68607),(0, 0.33116, 0.68178),(0, 0.33582, 0.67746),
(0, 0.34046, 0.6731),(0, 0.34509, 0.66871),(0, 0.3497, 0.66429),
(0, 0.3543, 0.65984),(0, 0.35888, 0.65536),(0, 0.36346, 0.65085),
(0, 0.36803, 0.6463),(0, 0.37258, 0.64173),(0, 0.37713, 0.63713),
(0, 0.38167, 0.6325),(0, 0.38618, 0.62783),(0, 0.39071, 0.62313),
(0, 0.39523, 0.6184),(0, 0.39972, 0.61363),(0, 0.40423, 0.60885),
(0, 0.40872, 0.60402),(0, 0.41321, 0.59915),(0, 0.41769, 0.59426),
(0, 0.42215, 0.58932),(0, 0.42663, 0.58437),(0, 0.4311, 0.57937),
(0, 0.43556, 0.57433),(0, 0.44001, 0.56926),(0, 0.44446, 0.56416),
(0, 0.44891, 0.55902),(0, 0.45334, 0.55384),(0, 0.45778, 0.54863),
(0, 0.46222, 0.54336),(0, 0.46665, 0.53802),(0, 0.47105, 0.53253),
(0, 0.47545, 0.52691),(0, 0.47982, 0.52115),(0, 0.48417, 0.51525),
(0, 0.48852, 0.50921),(0, 0.49284, 0.50301),(0, 0.49717, 0.49668),
(0, 0.50147, 0.49022),(0, 0.50575, 0.48359),(0, 0.51003, 0.47682),
(0, 0.51431, 0.4699),(0, 0.51856, 0.4628),(0.0097866, 0.52281, 0.45558),
(0.023896, 0.52704, 0.44818),(0.038383, 0.53126, 0.44061),(0.051763, 0.53547, 0.43289),
(0.063442, 0.53968, 0.42499),(0.073828, 0.54388, 0.41692),(0.083244, 0.54807, 0.40866),
(0.092062, 0.55225, 0.40022),(0.10019, 0.55642, 0.39159),(0.10786, 0.56059, 0.38276),
(0.11513, 0.56474, 0.37372),(0.12206, 0.56889, 0.36445),(0.12871, 0.57304, 0.35498),
(0.13507, 0.57718, 0.34524),(0.14115, 0.58131, 0.33527),(0.14697, 0.58544, 0.32499),
(0.15257, 0.58954, 0.31449),(0.15773, 0.59367, 0.30393),(0.16231, 0.59779, 0.29352),
(0.16631, 0.60191, 0.28332),(0.16984, 0.60603, 0.27332),(0.17292, 0.61015, 0.26352),
(0.17565, 0.61427, 0.25387),(0.17811, 0.6184, 0.24439),(0.18021, 0.62252, 0.23514),
(0.18207, 0.62664, 0.22606),(0.18374, 0.63076, 0.21715),(0.18522, 0.63487, 0.2084),
(0.18649, 0.63898, 0.19982),(0.18765, 0.64309, 0.19148),(0.18863, 0.6472, 0.18334),
(0.18951, 0.6513, 0.1754),(0.19029, 0.65539, 0.16767),(0.19099, 0.65948, 0.16013),
(0.19162, 0.66357, 0.15293),(0.19219, 0.66765, 0.14604),(0.19272, 0.67172, 0.13937),
(0.1932, 0.67579, 0.13311),(0.19366, 0.67985, 0.1272),(0.1941, 0.68391, 0.1218),
(0.19451, 0.68797, 0.11692),(0.1949, 0.69202, 0.11259),(0.19529, 0.69606, 0.10881),
(0.19569, 0.70009, 0.10581),(0.19611, 0.70412, 0.10358),(0.19653, 0.70814, 0.10209),
(0.19694, 0.71215, 0.10139),(0.19736, 0.71617, 0.10116),(0.19779, 0.72018, 0.10101),
(0.19823, 0.7242, 0.10087),(0.19868, 0.72822, 0.10073),(0.19914, 0.73225, 0.1006),
(0.19961, 0.73627, 0.10048),(0.20009, 0.74031, 0.10036),(0.20058, 0.74434, 0.10025),
(0.20108, 0.74838, 0.10015),(0.20159, 0.75242, 0.10006),(0.20211, 0.75647, 0.099977),
(0.20265, 0.76052, 0.099902),(0.2032, 0.76457, 0.099835),(0.20376, 0.76862, 0.099777),
(0.20433, 0.77267, 0.099729),(0.20488, 0.77674, 0.099691),(0.20546, 0.7808, 0.099663),
(0.20608, 0.78487, 0.099645),(0.20669, 0.78894, 0.099637),(0.20729, 0.79301, 0.099641),
(0.20791, 0.79708, 0.099656),(0.20855, 0.80116, 0.099683),(0.2092, 0.80523, 0.09972),
(0.20987, 0.80932, 0.09977),(0.21055, 0.8134, 0.099833),(0.21125, 0.81749, 0.099908),
(0.21196, 0.82159, 0.099996),(0.21266, 0.82567, 0.1001),(0.2134, 0.82977, 0.10021),
(0.21454, 0.83386, 0.10034),(0.21746, 0.83784, 0.10049),(0.22334, 0.84166, 0.10065),
(0.23182, 0.84531, 0.10083),(0.24228, 0.84884, 0.10102),(0.25428, 0.85224, 0.10123),
(0.26735, 0.85553, 0.10145),(0.28125, 0.85872, 0.10168),(0.29571, 0.86184, 0.10193),
(0.31067, 0.86485, 0.10219),(0.32594, 0.86779, 0.10246),(0.34137, 0.87065, 0.10275),
(0.35684, 0.87346, 0.10306),(0.3724, 0.87619, 0.10337),(0.38805, 0.87886, 0.10369),
(0.40366, 0.88147, 0.10401),(0.41921, 0.88402, 0.10433),(0.43478, 0.88651, 0.10468),
(0.45028, 0.88894, 0.10507),(0.4657, 0.89133, 0.10549),(0.48111, 0.89365, 0.10592),
(0.49641, 0.89593, 0.10637),(0.51169, 0.89815, 0.10681),(0.52691, 0.90032, 0.10725),
(0.54202, 0.90245, 0.1077),(0.55713, 0.90452, 0.10817),(0.57212, 0.90655, 0.10869),
(0.5871, 0.90852, 0.10927),(0.60201, 0.91045, 0.10983),(0.61689, 0.91232, 0.11037),
(0.63169, 0.91415, 0.11095),(0.64646, 0.91593, 0.11155),(0.66118, 0.91767, 0.11216),
(0.67585, 0.91934, 0.11281),(0.69049, 0.92098, 0.11347),(0.70508, 0.92257, 0.11409),
(0.71966, 0.92411, 0.11477),(0.73418, 0.9256, 0.11552),(0.74868, 0.92704, 0.11627),
(0.76314, 0.92844, 0.11699),(0.77759, 0.92979, 0.11777),(0.79201, 0.93109, 0.11859),
(0.8064, 0.93233, 0.11937),(0.82077, 0.93353, 0.12019),(0.83511, 0.93468, 0.12102),
(0.84946, 0.93578, 0.12189),(0.86375, 0.93684, 0.12278),(0.87808, 0.93783, 0.1237),
(0.89234, 0.93878, 0.12464),(0.90664, 0.93968, 0.12562),(0.92088, 0.94052, 0.12657),
(0.93514, 0.94131, 0.12755),(0.94939, 0.94206, 0.12857),(0.9636, 0.94275, 0.12961),
(0.97785, 0.94338, 0.13068),(0.99205, 0.94397, 0.13172),(1, 0.94449, 0.13281),
(1, 0.94497, 0.13392),(1, 0.94539, 0.13505),(1, 0.94574, 0.13614),(1, 0.94606, 0.13735),
(0.96837, 0.97701, 0.035491),(0.96656, 0.97219, 0.051675),(0.96476, 0.96736, 0.065035),
(0.96302, 0.96252, 0.076202),(0.96132, 0.95766, 0.086007),(0.95966, 0.95279, 0.094821),
(0.95807, 0.94791, 0.10259),(0.95649, 0.94302, 0.1098),(0.955, 0.93811, 0.11627),
(0.95354, 0.93319, 0.12222),(0.95212, 0.92826, 0.12782),(0.95079, 0.92331, 0.13295),
(0.94948, 0.91834, 0.13776),(0.94825, 0.91336, 0.14216),(0.94707, 0.90836, 0.14627),
(0.94594, 0.90334, 0.15015),(0.9449, 0.8983, 0.15365),(0.94391, 0.89324, 0.15697),
(0.94298, 0.88816, 0.16003),(0.94214, 0.88305, 0.16286),(0.94135, 0.87792, 0.16542),
(0.94064, 0.87277, 0.16783),(0.94002, 0.86758, 0.16998),(0.93948, 0.86238, 0.17189),
(0.93901, 0.85713, 0.17363),(0.93863, 0.85186, 0.17517),(0.93835, 0.84655, 0.17643),
(0.93816, 0.84121, 0.17756),(0.93806, 0.83583, 0.17845),(0.93807, 0.83042, 0.17911),
(0.9382, 0.82494, 0.17952),(0.93845, 0.81943, 0.17972),(0.93881, 0.81387, 0.1797),
(0.93929, 0.80825, 0.17943),(0.93991, 0.80258, 0.17893),(0.94067, 0.79684, 0.17816),
(0.94159, 0.79104, 0.17704),(0.94268, 0.78516, 0.17562),(0.94394, 0.7792, 0.1739),
(0.9454, 0.77315, 0.1718),(0.94705, 0.767, 0.16936),(0.94892, 0.76075, 0.16643),
(0.95104, 0.75439, 0.16306),(0.95333, 0.74791, 0.15928),(0.95563, 0.74139, 0.15542),
(0.95785, 0.73488, 0.15162),(0.96, 0.72836, 0.14791),(0.96208, 0.72182, 0.14425),
(0.96409, 0.71529, 0.14067),(0.96603, 0.70874, 0.13715),(0.96791, 0.70217, 0.1337),
(0.9697, 0.69561, 0.13037),(0.97144, 0.68904, 0.12703),(0.9731, 0.68244, 0.12383),
(0.9747, 0.67585, 0.12072),(0.97623, 0.66924, 0.11771),(0.97771, 0.66262, 0.11474),
(0.97912, 0.65599, 0.11195),(0.98046, 0.64934, 0.10919),(0.98175, 0.64268, 0.10656),
(0.98297, 0.63601, 0.10396),(0.98414, 0.62932, 0.1015),(0.98524, 0.62262, 0.099162),
(0.98629, 0.6159, 0.096892),(0.98728, 0.60918, 0.094848),(0.98821, 0.60242, 0.092755),
(0.98909, 0.59566, 0.090908),(0.98991, 0.58886, 0.089115),(0.99067, 0.58206, 0.087485),
(0.99138, 0.57524, 0.085856),(0.99204, 0.5684, 0.084524),(0.99265, 0.56153, 0.083158),
(0.99322, 0.55463, 0.081958),(0.99375, 0.54769, 0.08063),(0.99426, 0.54072, 0.079325),
(0.99474, 0.53369, 0.078024),(0.9952, 0.52664, 0.076731),(0.99563, 0.51953, 0.075448),
(0.99603, 0.51237, 0.074177),(0.9964, 0.50517, 0.072849),(0.99675, 0.49792, 0.071543),
(0.99707, 0.49061, 0.070164),(0.99736, 0.48323, 0.068891),(0.99763, 0.47581, 0.067489),
(0.99787, 0.46831, 0.066112),(0.99808, 0.46074, 0.064752),(0.99827, 0.45312, 0.063411),
(0.99844, 0.4454, 0.061996),(0.99857, 0.43763, 0.060625),(0.99869, 0.42975, 0.059296),
(0.99878, 0.42178, 0.057786),(0.99884, 0.41373, 0.056502),(0.99887, 0.40558, 0.055042),
(0.99889, 0.39731, 0.05362),(0.99888, 0.38894, 0.052219),(0.99884, 0.38043, 0.050843),
(0.99878, 0.37178, 0.049508),(0.99869, 0.36301, 0.048013),(0.99859, 0.35407, 0.04666),
(0.99844, 0.34498, 0.045261),(0.9982, 0.3358, 0.043828),(0.99777, 0.32662, 0.042522),
(0.99713, 0.31749, 0.041278),(0.99629, 0.3084, 0.039989),(0.99526, 0.29935, 0.038765),
(0.99402, 0.29032, 0.0376),(0.99259, 0.28137, 0.036475),(0.99095, 0.27247, 0.035441),
(0.98911, 0.26367, 0.034197),(0.98708, 0.25493, 0.03322),(0.98485, 0.24624, 0.032279),
(0.98243, 0.23769, 0.031373),(0.97982, 0.22918, 0.0305),(0.97702, 0.22083, 0.029658),
(0.97403, 0.21254, 0.028843),(0.97086, 0.20443, 0.028053),(0.9675, 0.19646, 0.027287),
(0.96394, 0.18866, 0.026541),(0.9602, 0.18101, 0.025814),(0.95628, 0.17361, 0.025105),
(0.95218, 0.16642, 0.024411),(0.94791, 0.15946, 0.023731),(0.94345, 0.15278, 0.023063),
(0.93881, 0.14642, 0.022407),(0.93399, 0.14044, 0.02176),(0.92898, 0.13485, 0.021122),
(0.92381, 0.12964, 0.020493),(0.91846, 0.12485, 0.019871),(0.91299, 0.12043, 0.019262),
(0.90744, 0.11624, 0.018679),(0.90186, 0.11205, 0.018127),(0.89625, 0.10787, 0.017605),
(0.89061, 0.10383, 0.017114),(0.88494, 0.099784, 0.016653),(0.87924, 0.095821, 0.016218),
(0.87351, 0.09196, 0.015811),(0.86776, 0.088098, 0.015429),(0.86199, 0.084349, 0.015072),
(0.85618, 0.080598, 0.014739),(0.85034, 0.076937, 0.014432),(0.84448, 0.073512, 0.014147),
(0.83859, 0.070023, 0.013884),(0.83267, 0.066671, 0.01364),(0.82673, 0.063454, 0.013414),
(0.82076, 0.060328, 0.013209),(0.81476, 0.057256, 0.013025),(0.80874, 0.054357, 0.012862),
(0.80269, 0.051607, 0.01272),(0.79662, 0.049121, 0.012597),(0.79052, 0.046616, 0.012491),
(0.7844, 0.044355, 0.0124),(0.77825, 0.042283, 0.012323),(0.77207, 0.040359, 0.012259),
(0.76588, 0.038639, 0.012207),(0.75965, 0.037152, 0.012168),(0.7534, 0.035909, 0.012143),
(0.74714, 0.034746, 0.01213),(0.74086, 0.033657, 0.012129),(0.73461, 0.032755, 0.012137),
(0.72835, 0.031827, 0.012155),(0.7221, 0.030928, 0.012182),(0.71588, 0.030041, 0.012218),
(0.70965, 0.02917, 0.012261),(0.70343, 0.028316, 0.012312),(0.69722, 0.027476, 0.012369),
(0.69103, 0.026652, 0.012431),(0.68485, 0.025842, 0.012477),(0.67867, 0.025046, 0.012503),
(0.67251, 0.024264, 0.012512),(0.66635, 0.023496, 0.012501),(0.66022, 0.022741, 0.012473),
(0.65409, 0.022, 0.012427),(0.64796, 0.021273, 0.012363),(0.64186, 0.020559, 0.01228),
(0.63576, 0.019858, 0.012179),(0.62967, 0.019169, 0.012055),(0.62361, 0.018494, 0.011906),
(0.61754, 0.017832, 0.011726),(0.61149, 0.017183, 0.011511),(0.60544, 0.016546, 0.011263),
(0.59942, 0.015921, 0.010993),(0.59341, 0.015309, 0.010717),(0.5874, 0.014709, 0.010451),
(0.5814, 0.014124, 0.010192),(0.57543, 0.013545, 0.009933),(0.56945, 0.012988, 0.0096743),
(0.56349, 0.012487, 0.0094284),(0.55753, 0.012012, 0.0092082),(0.55158, 0.011488, 0.0090195),
(0.54563, 0.010949, 0.0088628),(0.53967, 0.010498, 0.0087375),(0.53373, 0.010138, 0.0086438),
(0.5278, 0.0098238, 0.0085815),(0.52187, 0.0095265, 0.0085507),(0.51595, 0.0092456, 0.0085516),
(0.51004, 0.008987, 0.0085858),(0.50414, 0.0087542, 0.0086553),(0.49823, 0.0085472, 0.008759),
(0.49232, 0.0083637, 0.0088963),(0.48643, 0.0082023, 0.0090674),(0.48054, 0.0080616, 0.0092717),
(0.47467, 0.0079404, 0.0095068),(0.46879, 0.0078379, 0.0097685),(0.46291, 0.0077536, 0.010055),
(0.45705, 0.0076872, 0.010385),(0.4512, 0.0076373, 0.010796),(0.44534, 0.0076027, 0.011316),
(0.4395, 0.0075827, 0.011859),(0.43365, 0.0075764, 0.012347),(0.42783, 0.0075831, 0.012837),
(0.42199, 0.007602, 0.01339),(0.41618, 0.0076329, 0.013987),(0.41036, 0.0076915, 0.014615),
(0.40447, 0.0079151, 0.015288),(0.39845, 0.0085386, 0.016013),(0.39225, 0.0096476, 0.016788),
(0.38589, 0.01125, 0.017611),(0.37938, 0.01319, 0.01848),(0.37272, 0.015331, 0.01939),
(0.36595, 0.017697, 0.02034),(0.35905, 0.020248, 0.021325),(0.35205, 0.022943, 0.022345),
(0.34497, 0.025745, 0.023396),(0.33779, 0.028628, 0.024479),(0.33052, 0.031576, 0.025595),
(0.32318, 0.034637, 0.02674),(0.31576, 0.037722, 0.02791),(0.30827, 0.040658, 0.029105),
(0.30072, 0.043344, 0.030323),(0.29312, 0.046013, 0.031565),(0.28544, 0.048622, 0.032846),
(0.27773, 0.050917, 0.034117),(0.26998, 0.053082, 0.035663),(0.26215, 0.055158, 0.036964),
(0.25428, 0.057101, 0.038324),(0.24634, 0.058898, 0.039699),(0.23838, 0.060521, 0.04109),
(0.23036, 0.062012, 0.042429),(0.2223, 0.063428, 0.043725),(0.2142, 0.064648, 0.045151),
(0.20605, 0.065769, 0.046454),(0.1978, 0.066788, 0.047758),(0.18954, 0.067685, 0.049199),
(0.18118, 0.068484, 0.050432),(0.17278, 0.069179, 0.051712),(0.16433, 0.069694, 0.052994),
(0.15576, 0.07009, 0.054295),(0.14709, 0.070382, 0.055604),(0.13843, 0.07042, 0.05692),
(0.12971, 0.070185, 0.058154),(0.12092, 0.069706, 0.059585),(0.112, 0.068994, 0.060816),
(0.10228, 0.068214, 0.06216),(0.091664, 0.06747, 0.063475),(0.07965, 0.066744, 0.064745),
(0.066001, 0.066024, 0.066019)
]}}
