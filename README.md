# xTools
A Mathematica package for working with GR and AdS/CFT.

Containing subpackages:

## xTension
Including xAct related functions. Especially functions that makes it easier for working with metric decompositions.

## Some notes on xDecomp
*xCoba* is excellent for component computation in xAct, but it cannot handle the case where the coordinate is partially defined, e.g., the spherically symmetric metric in $D$ dimension
$$
    \mathrm ds^2 = -h(r)\mathrm d t^2 + \frac{1}{f(r)}\mathrm d r^2 + r^2\mathrm d\Sigma_{k, D - 2}^2
$$
where $\mathrm d\Sigma_{k, D - 2}$ is the metric of maximally symmetric space with curvature $k$.

To handle such problems we can decompose the manifold into two manifolds as $M = M_1 \otimes M_2$ and introduce coordinates on $M_1$ only. The delta tensor $\delta_a^b$ can then be decomposed into the projection tensors onto $M_1$ and $M_2$
$$
    \delta_a^b = \delta^{(1), b}_a + \delta^{(2), b}_a
$$
Using these tensors we can project any tensor in $TM$ into $TM_n$, or extend any tensor in $TM_n$ into $TM$. For convenience we'll adopt the following notations
$$
    \begin{gather*}
        T^{[n, a]\cdots}{}_\cdots \equiv T^{b\cdots}{}_\cdots \delta^{(n), a}_b, \, \omega^{(n)}_a \equiv \omega_{[n, a]}
    \end{gather*}
$$

We can define *generalized* ordinary derivative operator $\tilde\nabla_a$, satisfying
$$
    \tilde\nabla_a v^b = \tilde\nabla_a v^{(1), c}\delta^{(1), b}_c + \tilde\nabla_a v^{(2), c}\delta^{(2), b}_c
$$
Similarly for general tensors. We wish $\tilde\nabla^{(1)}_a$ and $\tilde\nabla^{(2)}_a$ can be regarded as derivative operators of $M_1$ and $M_2$ respectively, i.e., the Christoffel tensor relating a different operator $\tilde\nabla'^{(1)}_a$ and $\tilde\nabla^{(1)}_a$ involves only $M_1$ components. Firstly we note that only diagonal components of Christoffel tensor are involved
$$
    \tilde\nabla'^{(1)}_a v^b = \tilde\nabla^{(1)}_a v^b + \Gamma(\tilde\nabla', \tilde\nabla)_{ac}^d v^{(1), c} \delta^{(1), b}_d + \Gamma(\tilde\nabla', \tilde\nabla)_{ac}^d v^{(2), c} \delta^{(2), b}_d
$$
The third term can be eliminated by imposing
$$
    [\tilde\nabla'^{(1)}_a - \tilde\nabla^{(1)}_a, \tilde\nabla^{(2)}_b]f = 0
$$

Let $\tilde\partial_a$ to be the connection with vanishing curvature, again satisfying $[\tilde\nabla^{(1)}_a - \tilde\partial^{(1)}_a, \tilde\nabla^{(2)}_b]f = 0$, the Riemann tensor associated with $\tilde\nabla_a$ can also be computed as follows: First we decompose the torsion
$$
    \begin{aligned}
        [\tilde\nabla_a, \tilde\nabla_b] f & = -T(\tilde\nabla)_{ab}^c\tilde\nabla_c f\\
        & = [\tilde\nabla^{(1)}_a + \tilde\nabla^{(2)}_a, \tilde\nabla^{(1)}_b + \tilde\nabla^{(2)}_b]f\\
        & = -T(\tilde\nabla^{(1)})_{ab}^c\tilde\nabla_c f - T(\tilde\nabla^{(2)})_{ab}^c\tilde\nabla_c f + 2[\tilde\nabla^{(1)}_{[a}, \tilde\nabla^{(2)}_{b]}]f\\
        T(\tilde\nabla)_{ab}^c & = T(\tilde\nabla^{(1)})_{ab}^c + T(\tilde\nabla^{(2)})_{ab}^c + T(\tilde\nabla^{(1)}, \tilde\nabla^{(2)})_{ab}^c
    \end{aligned}
$$
so that $R(\tilde\nabla)_{abc}{}^d$ may be computed from
$$
    \begin{aligned}
        R(\tilde\nabla)_{abc}{}^d\omega_d & = [\tilde\nabla_a, \tilde\nabla_b]\omega_c + T(\tilde\nabla)_{ab}^d \tilde\nabla_d\omega_c\\
        & = [\tilde\nabla^{(1)}_a + \tilde\nabla^{(2)}_a, \tilde\nabla^{(1)}_b + \tilde\nabla^{(2)}_b](\omega^{(1)}_c + \omega^{(2)}_c)+ T(\tilde\nabla)_{ab}^d \tilde\nabla_d\omega_c\\
        & = [\tilde\nabla^{(1)}_a, \tilde\nabla^{(1)}_b]\omega^{(1)}_c + [\tilde\nabla^{(2)}_a, \tilde\nabla^{(2)}_b]\omega^{(2)}_c+ T(\tilde\nabla)_{ab}^d \tilde\nabla_d\omega_c\\
        & \quad - T(\tilde\nabla^{(1)})_{ab}^d\tilde\nabla^{(1)}_d\omega^{(2)}_c - T(\tilde\nabla^{(2)})_{ab}^d\tilde\nabla^{(2)}_d\omega^{(1)}_c + 2[\tilde\nabla^{(1)}_{[a}, \tilde\nabla^{(2)}_{b]}]\omega_c\\
        & = R(\tilde\nabla^{(1)})_{abc}{}^d\omega_d + R(\tilde\nabla^{(2)})_{abc}{}^d\omega_d + 2[\tilde\nabla^{(1)}_{[a}, \tilde\nabla^{(2)}_{b]}]\omega_c + T(\tilde\nabla^{(1)}, \tilde\nabla^{(2)})_{ab}^d\tilde\nabla_d\omega_c\\
        & \equiv R(\tilde\nabla^{(1)})_{abc}{}^d\omega_d + R(\tilde\nabla^{(2)})_{abc}{}^d\omega_d + R(\tilde\nabla^{(1)}, \tilde\nabla^{(2)})_{abc}{}^d\omega_d
    \end{aligned}
$$
work out the cross term, in terms of $\tilde\partial_a$
$$
    \begin{aligned}
        2[\tilde\partial^{(1)}_{[a}, \tilde\partial^{(2)}_{b]}]\omega_c & = -T(\tilde\partial^{(1)}, \tilde\partial^{(2)})_{ab}^d\tilde\partial_d\omega_c\\

        2[\tilde\nabla^{(1)}_{[a}, \tilde\nabla^{(2)}_{b]}]\omega^{(1)}_c & = 2\tilde\nabla^{(1)}_{[a}\tilde\nabla^{(2)}_{b]}\omega^{(1)}_c - 2\tilde\nabla^{(2)}_{[b}\tilde\nabla^{(1)}_{a]}\omega^{(1)}_c\\
        & = 2\tilde\nabla^{(1)}_{[a}\tilde\partial^{(2)}_{b]}\omega^{(1)}_c - 2\tilde\partial^{(2)}_{[b}\left[\tilde\partial^{(1)}_{a]}\omega^{(1)}_c - \Gamma(\tilde\nabla^{(1)}, \tilde\partial^{(1)})_{a]c}^d\omega^{(1)}_d\right]\\
        & = 2[\tilde\partial^{(1)}_{[a}, \tilde\partial^{(2)}_{b]}]\omega^{(1)}_c - 2\tilde\partial^{(2)}_{[a}\Gamma(\tilde\nabla^{(1)}, \tilde\partial^{(1)})_{b]c}^d\omega^{(1)}_d\\

        2[\tilde\nabla^{(1)}_{[a}, \tilde\nabla^{(2)}_{b]}]\omega_c & = 2[\tilde\partial^{(1)}_{[a}, \tilde\partial^{(2)}_{b]}]\omega_c - 2\tilde\partial^{(2)}_{[a}\Gamma(\tilde\nabla^{(1)}, \tilde\partial^{(1)})_{b]c}^d\omega^{(1)}_d - 2\tilde\partial^{(1)}_{[a}\Gamma(\tilde\nabla^{(2)}, \tilde\partial^{(2)})_{b]c}^d\omega^{(2)}_d\\

        R(\tilde\nabla^{(1)}, \tilde\nabla^{(2)})_{abc}{}^d\omega_ d & = - 2\tilde\partial^{(2)}_{[a}\Gamma(\tilde\nabla^{(1)}, \tilde\partial^{(1)})_{b]c}^d\omega^{(1)}_d - 2\tilde\partial^{(1)}_{[a}\Gamma(\tilde\nabla^{(2)}, \tilde\partial^{(2)})_{b]c}^d\omega^{(2)}_d + T(\tilde\nabla^{(1)}, \tilde\nabla^{(2)})_{ab}^d(\tilde\nabla_d - \tilde\partial_d)\omega_c
    \end{aligned}
$$
so that
$$
    R(\tilde\nabla^{(1)}, \tilde\nabla^{(2)})_{abc}{}^d = - 2\tilde\partial^{(2)}_{[a}\Gamma(\tilde\nabla^{(1)}, \tilde\partial^{(1)})_{b]c}^d - 2\tilde\partial^{(1)}_{[a}\Gamma(\tilde\nabla^{(2)}, \tilde\partial^{(2)})_{b]c}^d - T(\tilde\partial^{(1)}, \tilde\partial^{(2)})_{ab}^e\left[\Gamma(\tilde\nabla^{(1)}, \tilde\partial^{(1)})_{ec}^d + \Gamma(\tilde\nabla^{(2)}, \tilde\partial^{(2)})_{ec}^d\right]
$$

To compute the curvature of the real connection $\nabla_a$, we first compute the Christoffel tensor $\Gamma(\nabla, \tilde\nabla)_{ab}^c$, in the case of Levi-Civita connection, it is given by
$$
    \Gamma(\nabla, \tilde\nabla)_{ab}^c = \frac12 g^{cd}\left(\tilde\nabla_a g_{db} + \tilde\nabla_b g_{ad} - \tilde\nabla_d g_{ab}\right)
$$

The Riemann tensor is then
$$
    R(\nabla)_{abc}{}^d = R(\tilde\nabla)_{abc}{}^d -2\tilde\nabla_{[a} \Gamma(\nabla, \tilde\nabla)_{b]c}^d + 2\Gamma(\nabla, \tilde\nabla)_{c[a}^e\Gamma(\nabla, \tilde\nabla)_{b]e}^d
$$