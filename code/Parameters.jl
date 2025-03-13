# Parameters for calculating GCD, R0, N0
const (KJ, rhoJ, muJ, varPhi) = [0.75 * 1E3, (1/(12 * 1440)), (1 / (20 * 1440)), (3 / 1440)]
const (gV, gR, mu) = [(1/(5 * 1440)), (1/(2 * 1440)), (1/(21 * 1440))]
const (muH, KH) = [1/(365.25 * 65 * 1440), 1E3]
# const pQ = 1.0f0
const (bH, bB, eta, gH) = [0.5f0, 0.5f0, 1 / (7 * 1440), 1 / (5 * 1440)]