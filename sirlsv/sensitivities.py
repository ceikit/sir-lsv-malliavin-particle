"""The four auxiliary processes of the paper (integrated Malliavin sensitivities).

With g = S dsigma/dS, v_i the Hull-White variances and W the Q^d Brownian motion of the spot:

    dl_d = l_d [ g (dW - sigma dt) - k_d dt ] + [ v_d - z_f + rho_Sd s_d sigma ] dt
    dz_f = [ -(k_d + k_f) z_f - rho_Sf s_f g l_d + rho_df s_d s_f ] dt
    dl_f = l_f [ g dW - k_f dt ]               + [ z_d - v_f + rho_Sf s_f sigma ] dt
    dz_d = [ -(k_d + k_f) z_d + rho_Sd s_d g l_f + rho_df s_d s_f ] dt

all starting at 0. The calibrated local variance is
sigma_Dupire^2 - 2 E[ l_d - l_f | S_t = K ] under the domestic t-forward measure.
"""


def sensitivity_step(ld, zf, lf, zd, g, sig, dW, dt, vd, vf, p):
    """One Euler step of (l_d, z_f, l_f, z_d).

    g, sig : S dsigma/dS and sigma at the start of the step, per particle
    dW     : increment of the Q^d Brownian motion of the spot, per particle
    vd, vf : Hull-White variances v_d(t), v_f(t) at the start of the step
    p      : RateParams
    """
    kd, sd, kf, sf, rSd, rSf, rdf = p.kd, p.sd, p.kf, p.sf, p.rSd, p.rSf, p.rdf
    ld_n = ld * (1 + g * (dW - sig * dt) - kd * dt) + (vd - zf + rSd * sd * sig) * dt
    zf_n = zf + (-(kd + kf) * zf - rSf * sf * g * ld + rdf * sd * sf) * dt
    lf_n = lf * (1 + g * dW - kf * dt) + (zd - vf + rSf * sf * sig) * dt
    zd_n = zd + (-(kd + kf) * zd + rSd * sd * g * lf + rdf * sd * sf) * dt
    return ld_n, zf_n, lf_n, zd_n
