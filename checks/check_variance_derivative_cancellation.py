"""Symbolic check of the two expanded ell PDEs; not an existence test."""
import sympy as sp

x, u, v = sp.symbols('x u v', real=True)
kd, kf, beta, gamma, c, nd, nf, fd, ff = sp.symbols(
    'kd kf beta gamma c nd nf fd ff', real=True)
a = sp.Function('a', positive=True)(x)
m = sp.Function('m')(x, u, v)
s = sp.sqrt(a)
g = sp.diff(s, x)
rd, rf = fd + u, ff + v
b = [rd-rf-a/2, -kd*u, -kf*v-gamma*s]
z = [x, u, v]
q = sp.Matrix([[a, beta*s, gamma*s], [beta*s, nd**2, c],
               [gamma*s, c, nf**2]])
killed = -sum(sp.diff(bi*m, zi) for bi, zi in zip(b, z)) - rd*m
killed += sum(sp.diff(q[i, j]*m, z[i], z[j])/2
              for i in range(3) for j in range(3))
flux = -sum(sp.diff(g*hi*m, zi) for hi, zi in zip([s, beta, gamma], z))
principal = sum(q[i, j]*sp.diff(m, z[i], z[j])/2
                for i in range(3) for j in range(3))
expanded = principal + (-rd+rf+a/2+sp.diff(a, x)/2)*sp.diff(m, x)
expanded += kd*u*sp.diff(m, u)+(kf*v+gamma*s)*sp.diff(m, v)
assert sp.simplify(killed+flux-(kd+s*g)*m
                   -expanded-(kf-rd)*m) == 0
assert sp.simplify(killed+flux-kf*m
                   -expanded-(kd-rd+sp.diff(a, x)/2)*m) == 0
second = sp.Symbol('second')
passive = sp.expand(killed).subs(sp.diff(a, x, 2), second)
assert sp.simplify(sp.diff(passive, second)-m/2) == 0
print('PASS: domestic and foreign ell cancellations; passive a_xx coefficient = m/2')
