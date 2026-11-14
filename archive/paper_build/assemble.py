"""Assemble the merged paper from the original 2016 draft (Sections 1-2, the proof of the
calibration condition, the Malliavin primer) and the rewritten Section 3.

usage: python assemble.py <original_draft.txt> <output.tex> <check_sensitivities_output.txt>
The original draft and the results log are not part of this repository.
"""
import re
import sys

src_path, out_path, mc_path = sys.argv[1], sys.argv[2], sys.argv[3]
src = open(src_path, encoding="utf-8").read().split("\n")
L = lambda a, b: "\n".join(src[a - 1:b])          # 1-indexed inclusive line range
import os
HERE = os.path.dirname(os.path.abspath(__file__))
here = lambda name: open(os.path.join(HERE, name), encoding="utf-8").read()


def sub1(text, old, new):
    assert text.count(old) == 1, (text.count(old), old)
    return text.replace(old, new)


# ---------------------------------------------------------------- preamble
pre = L(3, 46)
pre = sub1(pre, r"\usepackage{graphicx}", "\\usepackage{graphicx}\n\\usepackage{booktabs}")
pre += "\n\\newcommand{\\MD}{\\text{\\DH}}\n\\date{Draft of 2016, revised November 2026}\n"

abstract = r"""
\begin{abstract}
We derive the particle (McKean) calibration of local volatility models with stochastic interest rates, following Guyon and Henry-Labord\`ere, and extend it to long-dated FX with Hull-White domestic and foreign rates. The three-factor calibration condition is written under a single measure, so that one system of particles suffices. Its Malliavin representation is obtained without solving the coupled two-parameter system of Malliavin derivatives: integrating the derivatives against the Hull-White kernels yields four auxiliary processes per particle, simulated forward together with the state, and a correction to the Dupire local variance which involves neither $\partial_K^2C$ nor tail expectations. The representation reduces to the known one-factor formula, matches the Gaussian closed form, and its building blocks are verified by Monte Carlo. A full calibration to a market smile has not yet been run, and existence of the calibrated model is not addressed.
\end{abstract}
"""

# ---------------------------------------------------------------- sections 1-2
s12 = L(54, 299)
s12 = sub1(s12, r"r_t^0 = \mathbb{E}^{\mathbb{Q}} [r_t] = -\partial_t\log P_{0t}",
           r"r_t^0 = \mathbb{E}^{\mathbb{Q}^t} [r_t] = -\partial_t\log P_{0t}")
s12 = sub1(s12, r"(1-e^{=2\kappa t})", r"(1-e^{-2\kappa t})")
s12 = sub1(s12, r"=B_s+\sigma_r\int_0^tA(s,T)ds", r"=B_t+\sigma_r\int_0^tA(s,T)ds")
s12 = sub1(s12, r"\Xi_t &= \sigma_t\frac{e^{2\kappa t}-1}{2\kappa}", r"\Xi_t &= \sigma_r\frac{e^{2\kappa t}-1}{2\kappa}")
s12 = sub1(s12, r"where $B^i$ and $W^i$ are $\mathbb{Q}$-Brownian motions.",
           r"where $B^i$ and $W^i$ are $\mathbb{Q}^T$-Brownian motions.")
s12 = sub1(s12, r"from (\ref{eq: hull_white_integral}):", r"from (\ref{eq: hullWhiteSolution}):")
s12 = sub1(s12, r"(W_{t_{j+1}}^T-W^T_{t_{j}})\int_{t_j}^{t_{j+1}} \frac{1-e^{-\kappa(T-t)}}{\kappa}dB_t^T",
           r"(W_{t_{j+1}}^T-W^T_{t_{j}})+\sigma_r\int_{t_j}^{t_{j+1}} \frac{1-e^{-\kappa(T-t)}}{\kappa}dB_t^T")
s12 = sub1(s12, r"\centerline{\includegraphics[scale=0.5]{volatility.jpg} }",
           r"\IfFileExists{volatility.jpg}{\centerline{\includegraphics[scale=0.5]{volatility.jpg}}}"
           r"{\centerline{\fbox{\parbox{0.6\textwidth}{\centering\small [figure \texttt{volatility.jpg}: place the file next to this source]}}}}")
s12 = sub1(s12, r"\item[1.] Initialize $j=1$ and set $\sigma_N(t,K)=\frac{\sigma_{Dup}(0,K)}{grid\ range}$",
           r"\item[1.] Initialize $j=1$ and set $\sigma_N(t,K)=\sigma_{Dup}(0,K)$ on the first time step")
s12 = sub1(s12, r"\subsection{Malliavin representation of the local volatility}",
           r"\subsection{Malliavin representation of the local volatility}\label{subsec: malliavinGeneral}")
s12 = sub1(s12, r"\subsection{Malliavin representation under Hull-White}",
           r"\subsection{Malliavin representation under Hull-White}\label{subsec: malliavin1F}")

proof = L(346, 397)
primer = sub1(L(649, 696), "three-factor FX SIR-SV model", "three-factor FX SIR-LV model")

# ---------------------------------------------------------------- Monte Carlo table
rows = {}
block = None
for line in open(mc_path):
    if "(I)" in line:
        block = "d"
    elif "(II)" in line:
        block = "f"
    m = re.match(r"\s+(\d\.\d\d) \| (\S+) (\S+)\s+(\S+) \| (\S+)\s+(\S+)", line)
    if m and block:
        K, direct, sens, z, naive, zn = m.groups()
        rows.setdefault(K, {})[block] = (float(direct) * 1e3, float(sens) * 1e3, float(z), float(naive) * 1e3, float(zn))
tab = []
fz = lambda z: "0.00" if abs(z) < 0.005 else f"{z:+.2f}"
for K in sorted(rows):
    d, f = rows[K]["d"], rows[K]["f"]
    tab.append(f"{K} & ${d[0]:.4f}$ & ${d[1]:.4f}$ & ${fz(d[2])}$ & ${f[0]:.4f}$ & ${f[1]:.4f}$ & ${fz(f[2])}$ & ${f[3]:.4f}$ & ${fz(f[4])}$\\\\")
assert len(tab) == 6
s3d = here("sec3d.tex")
start = s3d.index("0.70 & $2.425$")
end = s3d.index("\\bottomrule")
s3d = s3d[:start] + "\n".join(tab) + "\n" + s3d[end:]

body = "\n".join([
    pre, r"\begin{document}", "", r"\maketitle", abstract, r"\tableofcontents", "",
    s12, here("sec2_extra.tex"), here("sec3a.tex"), proof, here("sec3b.tex"), primer, here("sec3c.tex"), s3d])
body = re.sub(r"\\DH(?![a-zA-Z])", r"\\MD", body)
body = body.replace("\\newcommand{\\MD}{\\text{\\MD}}", "\\newcommand{\\MD}{\\text{\\DH}}")
open(out_path, "w", encoding="utf-8").write(body)
print("written", out_path, len(body.split("\n")), "lines")
