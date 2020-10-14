// Copyright 2019 Global Phasing Ltd.

// X-ray scattering factor coefficients from International Tables
// for Crystallography Volume C, edition from 1992 or later.
// The same data is in CCP4 (file atomsf.lib) and in cctbx (table="it1992").
// But here we have only one entry per element. The original table also
// lists ions, for example Fe2+ and Fe3+ in addition to Fe.
// Scattering factors from ITC are approximated by four Gaussians
// and a constant value - together 9 coefficients known as Cromer-Mann
// coefficients. The approximation holds for sin(theta)/lambda < 2 A^-1.

#ifndef GEMMI_IT92_HPP_
#define GEMMI_IT92_HPP_

#include "formfact.hpp"  // for GaussianCoef
#include "elem.hpp"      // for El

namespace gemmi {

#if defined(__GNUC__) && __GNUC__-0 > 4
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#endif

template<class Real>
struct IT92 {
  using Coef = GaussianCoef<4, 1, Real>;
  static const Coef data[98];

  static bool has(El el) {
    return el != El::X && (el <= El::Cf || el == El::D);
  }

  static const Coef& get(El el) {
    if (el == El::D)
      el = El::H;
    return data[(int)el - 1];
  }

  static const Coef* get_ptr(El el) {
    return has(el) ? &get(el) : nullptr;
  }
};

template<class Real>
const typename IT92<Real>::Coef IT92<Real>::data[98] = {
  // a1, a2, a3, a4, b1, b2, b3, b4, c
  {0.493002, 0.322912, 0.140191, 0.04081, 10.5109, 26.1257, 3.14236, 57.7997, 0.003038}, // H
  {0.8734, 0.6309, 0.3112, 0.178, 9.1037, 3.3568, 22.9276, 0.9821, 0.0064}, // He
  {1.1282, 0.7508, 0.6175, 0.4653, 3.9546, 1.0524, 85.3905, 168.261, 0.0377}, // Li
  {1.5919, 1.1278, 0.5391, 0.7029, 43.6427, 1.8623, 103.483, 0.542, 0.0385}, // Be
  {2.0545, 1.3326, 1.0979, 0.7068, 23.2185, 1.021, 60.3498, 0.1403, -0.1932}, // B
  {2.31, 1.02, 1.5886, 0.865, 20.8439, 10.2075, 0.5687, 51.6512, 0.2156}, // C
  {12.2126, 3.1322, 2.0125, 1.1663, 0.0057, 9.8933, 28.9975, 0.5826, -11.529}, // N
  {3.0485, 2.2868, 1.5463, 0.867, 13.2771, 5.7011, 0.3239, 32.9089, 0.2508}, // O
  {3.5392, 2.6412, 1.517, 1.0243, 10.2825, 4.2944, 0.2615, 26.1476, 0.2776}, // F
  {3.9553, 3.1125, 1.4546, 1.1251, 8.4042, 3.4262, 0.2306, 21.7184, 0.3515}, // Ne
  {4.7626, 3.1736, 1.2674, 1.1128, 3.285, 8.8422, 0.3136, 129.424, 0.676}, // Na
  {5.4204, 2.1735, 1.2269, 2.3073, 2.8275, 79.2611, 0.3808, 7.1937, 0.8584}, // Mg
  {6.4202, 1.9002, 1.5936, 1.9646, 3.0387, 0.7426, 31.5472, 85.0886, 1.1151}, // Al
  {6.2915, 3.0353, 1.9891, 1.541, 2.4386, 32.3337, 0.6785, 81.6937, 1.1407}, // Si
  {6.4345, 4.1791, 1.78, 1.4908, 1.9067, 27.157, 0.526, 68.1645, 1.1149}, // P
  {6.9053, 5.2034, 1.4379, 1.5863, 1.4679, 22.2151, 0.2536, 56.172, 0.8669}, // S
  {11.4604, 7.1964, 6.2556, 1.6455, 0.0104, 1.1662, 18.5194, 47.7784, -9.5574}, // Cl
  {7.4845, 6.7723, 0.6539, 1.6442, 0.9072, 14.8407, 43.8983, 33.3929, 1.4445}, // Ar
  {8.2186, 7.4398, 1.0519, 0.8659, 12.7949, 0.7748, 213.187, 41.6841, 1.4228}, // K
  {8.6266, 7.3873, 1.5899, 1.0211, 10.4421, 0.6599, 85.7484, 178.437, 1.3751}, // Ca
  {9.189, 7.3679, 1.6409, 1.468, 9.0213, 0.5729, 136.108, 51.3531, 1.3329}, // Sc
  {9.7595, 7.3558, 1.6991, 1.9021, 7.8508, 0.5, 35.6338, 116.105, 1.2807}, // Ti
  {10.2971, 7.3511, 2.0703, 2.0571, 6.8657, 0.4385, 26.8938, 102.478, 1.2199}, // V
  {10.6406, 7.3537, 3.324, 1.4922, 6.1038, 0.392, 20.2626, 98.7399, 1.1832}, // Cr
  {11.2819, 7.3573, 3.0193, 2.2441, 5.3409, 0.3432, 17.8674, 83.7543, 1.0896}, // Mn
  {11.7695, 7.3573, 3.5222, 2.3045, 4.7611, 0.3072, 15.3535, 76.8805, 1.0369}, // Fe
  {12.2841, 7.3409, 4.0034, 2.3488, 4.2791, 0.2784, 13.5359, 71.1692, 1.0118}, // Co
  {12.8376, 7.292, 4.4438, 2.38, 3.8785, 0.2565, 12.1763, 66.3421, 1.0341}, // Ni
  {13.338, 7.1676, 5.6158, 1.6735, 3.5828, 0.247, 11.3966, 64.8126, 1.191}, // Cu
  {14.0743, 7.0318, 5.1652, 2.41, 3.2655, 0.2333, 10.3163, 58.7097, 1.3041}, // Zn
  {15.2354, 6.7006, 4.3591, 2.9623, 3.0669, 0.2412, 10.7805, 61.4135, 1.7189}, // Ga
  {16.0816, 6.3747, 3.7068, 3.683, 2.8509, 0.2516, 11.4468, 54.7625, 2.1313}, // Ge
  {16.6723, 6.0701, 3.4313, 4.2779, 2.6345, 0.2647, 12.9479, 47.7972, 2.531}, // As
  {17.0006, 5.8196, 3.9731, 4.3543, 2.4098, 0.2726, 15.2372, 43.8163, 2.8409}, // Se
  {17.1789, 5.2358, 5.6377, 3.9851, 2.1723, 16.5796, 0.2609, 41.4328, 2.9557}, // Br
  {17.3555, 6.7286, 5.5493, 3.5375, 1.9384, 16.5623, 0.2261, 39.3972, 2.825}, // Kr
  {17.1784, 9.6435, 5.1399, 1.5292, 1.7888, 17.3151, 0.2748, 164.934, 3.4873}, // Rb
  {17.5663, 9.8184, 5.422, 2.6694, 1.5564, 14.0988, 0.1664, 132.376, 2.5064}, // Sr
  {17.776, 10.2946, 5.72629, 3.26588, 1.4029, 12.8006, 0.125599, 104.354, 1.91213}, // Y
  {17.8765, 10.948, 5.41732, 3.65721, 1.27618, 11.916, 0.117622, 87.6627, 2.06929}, // Zr
  {17.6142, 12.0144, 4.04183, 3.53346, 1.18865, 11.766, 0.204785, 69.7957, 3.75591}, // Nb
  {3.7025, 17.2356, 12.8876, 3.7429, 0.2772, 1.0958, 11.004, 61.6584, 4.3875}, // Mo
  {19.1301, 11.0948, 4.64901, 2.71263, 0.864132, 8.14487, 21.5707, 86.8472, 5.40428}, // Tc
  {19.2674, 12.9182, 4.86337, 1.56756, 0.80852, 8.43467, 24.7997, 94.2928, 5.37874}, // Ru
  {19.2957, 14.3501, 4.73425, 1.28918, 0.751536, 8.21758, 25.8749, 98.6062, 5.328}, // Rh
  {19.3319, 15.5017, 5.29537, 0.605844, 0.698655, 7.98929, 25.2052, 76.8986, 5.26593}, // Pd
  {19.2808, 16.6885, 4.8045, 1.0463, 0.6446, 7.4726, 24.6605, 99.8156, 5.179}, // Ag
  {19.2214, 17.6444, 4.461, 1.6029, 0.5946, 6.9089, 24.7008, 87.4825, 5.0694}, // Cd
  {19.1624, 18.5596, 4.2948, 2.0396, 0.5476, 6.3776, 25.8499, 92.8029, 4.9391}, // In
  {19.1889, 19.1005, 4.4585, 2.4663, 5.8303, 0.5031, 26.8909, 83.9571, 4.7821}, // Sn
  {19.6418, 19.0455, 5.0371, 2.6827, 5.3034, 0.4607, 27.9074, 75.2825, 4.5909}, // Sb
  {19.9644, 19.0138, 6.14487, 2.5239, 4.81742, 0.420885, 28.5284, 70.8403, 4.352}, // Te
  {20.1472, 18.9949, 7.5138, 2.2735, 4.347, 0.3814, 27.766, 66.8776, 4.0712}, // I
  {20.2933, 19.0298, 8.9767, 1.99, 3.9282, 0.344, 26.4659, 64.2658, 3.7118}, // Xe
  {20.3892, 19.1062, 10.662, 1.4953, 3.569, 0.3107, 24.3879, 213.904, 3.3352}, // Cs
  {20.3361, 19.297, 10.888, 2.6959, 3.216, 0.2756, 20.2073, 167.202, 2.7731}, // Ba
  {20.578, 19.599, 11.3727, 3.28719, 2.94817, 0.244475, 18.7726, 133.124, 2.14678}, // La
  {21.1671, 19.7695, 11.8513, 3.33049, 2.81219, 0.226836, 17.6083, 127.113, 1.86264}, // Ce
  {22.044, 19.6697, 12.3856, 2.82428, 2.77393, 0.222087, 16.7669, 143.644, 2.0583}, // Pr
  {22.6845, 19.6847, 12.774, 2.85137, 2.66248, 0.210628, 15.885, 137.903, 1.98486}, // Nd
  {23.3405, 19.6095, 13.1235, 2.87516, 2.5627, 0.202088, 15.1009, 132.721, 2.02876}, // Pm
  {24.0042, 19.4258, 13.4396, 2.89604, 2.47274, 0.196451, 14.3996, 128.007, 2.20963}, // Sm
  {24.6274, 19.0886, 13.7603, 2.9227, 2.3879, 0.1942, 13.7546, 123.174, 2.5745}, // Eu
  {25.0709, 19.0798, 13.8518, 3.54545, 2.25341, 0.181951, 12.9331, 101.398, 2.4196}, // Gd
  {25.8976, 18.2185, 14.3167, 2.95354, 2.24256, 0.196143, 12.6648, 115.362, 3.58324}, // Tb
  {26.507, 17.6383, 14.5596, 2.96577, 2.1802, 0.202172, 12.1899, 111.874, 4.29728}, // Dy
  {26.9049, 17.294, 14.5583, 3.63837, 2.07051, 0.19794, 11.4407, 92.6566, 4.56796}, // Ho
  {27.6563, 16.4285, 14.9779, 2.98233, 2.07356, 0.223545, 11.3604, 105.703, 5.92046}, // Er
  {28.1819, 15.8851, 15.1542, 2.98706, 2.02859, 0.238849, 10.9975, 102.961, 6.75621}, // Tm
  {28.6641, 15.4345, 15.3087, 2.98963, 1.9889, 0.257119, 10.6647, 100.417, 7.56672}, // Yb
  {28.9476, 15.2208, 15.1, 3.71601, 1.90182, 9.98519, 0.261033, 84.3298, 7.97628}, // Lu
  {29.144, 15.1726, 14.7586, 4.30013, 1.83262, 9.5999, 0.275116, 72.029, 8.58154}, // Hf
  {29.2024, 15.2293, 14.5135, 4.76492, 1.77333, 9.37046, 0.295977, 63.3644, 9.24354}, // Ta
  {29.0818, 15.43, 14.4327, 5.11982, 1.72029, 9.2259, 0.321703, 57.056, 9.8875}, // W
  {28.7621, 15.7189, 14.5564, 5.44174, 1.67191, 9.09227, 0.3505, 52.0861, 10.472}, // Re
  {28.1894, 16.155, 14.9305, 5.67589, 1.62903, 8.97948, 0.382661, 48.1647, 11.0005}, // Os
  {27.3049, 16.7296, 15.6115, 5.83377, 1.59279, 8.86553, 0.417916, 45.0011, 11.4722}, // Ir
  {27.0059, 17.7639, 15.7131, 5.7837, 1.51293, 8.81174, 0.424593, 38.6103, 11.6883}, // Pt
  {16.8819, 18.5913, 25.5582, 5.86, 0.4611, 8.6216, 1.4826, 36.3956, 12.0658}, // Au
  {20.6809, 19.0417, 21.6575, 5.9676, 0.545, 8.4484, 1.5729, 38.3246, 12.6089}, // Hg
  {27.5446, 19.1584, 15.538, 5.52593, 0.65515, 8.70751, 1.96347, 45.8149, 13.1746}, // Tl
  {31.0617, 13.0637, 18.442, 5.9696, 0.6902, 2.3576, 8.618, 47.2579, 13.4118}, // Pb
  {33.3689, 12.951, 16.5877, 6.4692, 0.704, 2.9238, 8.7937, 48.0093, 13.5782}, // Bi
  {34.6726, 15.4733, 13.1138, 7.02588, 0.700999, 3.55078, 9.55642, 47.0045, 13.677}, // Po
  {35.3163, 19.0211, 9.49887, 7.42518, 0.68587, 3.97458, 11.3824, 45.4715, 13.7108}, // At
  {35.5631, 21.2816, 8.0037, 7.4433, 0.6631, 4.0691, 14.0422, 44.2473, 13.6905}, // Rn
  {35.9299, 23.0547, 12.1439, 2.11253, 0.646453, 4.17619, 23.1052, 150.645, 13.7247}, // Fr
  {35.763, 22.9064, 12.4739, 3.21097, 0.616341, 3.87135, 19.9887, 142.325, 13.6211}, // Ra
  {35.6597, 23.1032, 12.5977, 4.08655, 0.589092, 3.65155, 18.599, 117.02, 13.5266}, // Ac
  {35.5645, 23.4219, 12.7473, 4.80703, 0.563359, 3.46204, 17.8309, 99.1722, 13.4314}, // Th
  {35.8847, 23.2948, 14.1891, 4.17287, 0.547751, 3.41519, 16.9235, 105.251, 13.4287}, // Pa
  {36.0228, 23.4128, 14.9491, 4.188, 0.5293, 3.3253, 16.0927, 100.613, 13.3966}, // U
  {36.1874, 23.5964, 15.6402, 4.1855, 0.511929, 3.25396, 15.3622, 97.4908, 13.3573}, // Np
  {36.5254, 23.8083, 16.7707, 3.47947, 0.499384, 3.26371, 14.9455, 105.98, 13.3812}, // Pu
  {36.6706, 24.0992, 17.3415, 3.49331, 0.483629, 3.20647, 14.3136, 102.273, 13.3592}, // Am
  {36.6488, 24.4096, 17.399, 4.21665, 0.465154, 3.08997, 13.4346, 88.4834, 13.2887}, // Cm
  {36.7881, 24.7736, 17.8919, 4.23284, 0.451018, 3.04619, 12.8946, 86.003, 13.2754}, // Bk
  {36.9185, 25.1995, 18.3317, 4.24391, 0.437533, 3.00775, 12.4044, 83.7881, 13.2674}, // Cf
};

#if defined(__GNUC__) && __GNUC__-0 > 4
#pragma GCC diagnostic pop
#endif

} // namespace gemmi
#endif
