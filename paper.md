---
title: "A first-order Madelung tie-break on the round four-sphere: Kelvin parity, convex order and the Thomas-Fermi field"
author: "Michaël Perreault"
date: "2026-11-25"
format: "article"
toc: false
header_footer_policy: "all"
footer: "© 2026 Michaël Perreault · CC BY 4.0"
pageof: true
date_footer: "YYYY-MM-DD"
header_title: "A first-order Madelung tie-break on the round four-sphere"
---

**Abstract.** In the Demkov-Ostrovsky model of the Madelung rule, the one-electron levels with equal $N=n+l$ are degenerate at zero energy, and our earlier paper [P1] showed that the model is the circle-invariant conformal Laplacian of the round four-sphere. The rule's second part, lower $n$ first within a group, is not explained by the model. We prove that within every group the laws of $|\cos\theta|$ ($\theta$ the polar angle) for consecutive $l$ are in strict increasing convex order, with an explicit closed form from Gegenbauer polynomials and the Laguerre inequality (Theorems 1 and 2). It is not first-order dominance: the distribution functions cross. Hence a perturbation of the Tietz field that depends on $r$ alone, with Kelvin-even weakening increasing and convex in $|\cos\theta|$, splits every group lower $n$ first at first order (Theorem 3, valid near the Tietz field). For the Thomas-Fermi field, whose weakening is not convex at the poles, Theorem 3 is not applied: a direct interval computation of the first-order slopes certifies the splitting in every group with $N\le17$, relative to its nearest Tietz field (one fitted scale, robust by estimate for $0.40\le a\le0.65$ at $N\le8$; Corollary 4). A numerical shooting check finds no level crossing inside a group along the whole path from the Tietz to the Thomas-Fermi field for $N\le12$ (Proposition 5); at the Thomas-Fermi end the groups 11 and 12 already overlap, so only the order inside a group is claimed. We do not derive real atoms, the spin factor, or the order between groups.

**Keywords.** Madelung rule; Demkov-Ostrovsky model; Gegenbauer polynomials; convex order; Thomas-Fermi equation; interval arithmetic.

```{=latex}
\begingroup\small
This work is licensed under the Creative Commons Attribution 4.0 International License (CC BY 4.0), https://creativecommons.org/licenses/by/4.0/. The accompanying code is MIT-licensed.
\par\endgroup
```

# Introduction {#p5.intro}

The Madelung rule has two parts: levels fill in order of increasing $N=n+l$, and levels with equal $N$ fill in order of increasing $n$ (decreasing $l$). Demkov and Ostrovsky modelled the first part by the zero-energy coupling problem for the Tietz weight [DO71, DO72]; its coupling constants depend on $N$ only, so every group is degenerate. Paper 1 [P1] (https://doi.org/10.5281/zenodo.24522163) identifies this problem with the circle-invariant conformal Laplacian of the round four-sphere, and the Kelvin inversion with the swap of its poles. That explains the degeneracy of a group. It leaves the order inside a group as an assumption of the model: nothing in the sphere picture says why $3p$ should lie below $4s$. This paper asks why lower $n$ comes first.

The idea is that of Goudsmit and Richards [GR]: a field slightly different from the Tietz field lifts the degeneracy, and the lift depends on $l$. We make this precise at first order. A perturbation $w\mapsto w(1+\epsilon g)$ with $g$ depending on $r$ alone shifts the level of the phase-free block $l$ in group $N$, to first order, by $\epsilon Z_N$ times the mean of the Kelvin-even weakening $h$ under the density $p_{2l}$ of $\cos\theta$ (Section 5). So the question is an order between the laws of $|\cos\theta|$ for consecutive $l$.

Roadmap. Section 3 shows why the first-order stochastic order fails, Section 4 proves the convex order (Theorems 1 and 2), Section 5 the first-order tie-break (Theorem 3), and Section 6 treats the Thomas-Fermi field (Corollary 4, certified for every group with $N\le17$, and the numerical Proposition 5 along the path from the Tietz field); it ends with a description of where the Thomas-Fermi field departs from the round sphere. Section 7 states the limits, and Section 8 lists the code.

In words: within this one-electron model, a field weaker than the Tietz field near the poles of the four-sphere (mainly the far pole, folded onto the nucleus end by Kelvin parity) raises the coupling of low-$l$ states more, because their $|\cos\theta|$ law lies further out in the increasing convex order; the high-$l$ state of a group therefore binds first, which is lower $n$ first.

We do not claim to derive real atoms, the spin factor, or the order between groups. The window $\epsilon_N$ of Theorem 3 depends on $N$ and is not bounded below; first order alone does not reach the Thomas-Fermi field, which is why Proposition 5 is a separate, numerical statement.

Credit. Goudsmit and Richards argue part B of the rule by a gradual transition between the Thomas-Fermi limit and the hydrogenic limit; that is the qualitative form of our statement, without perturbation theory or a theorem. Negadi and Kibler [NK92] (a $q$-deformation) reproduce the order with a fitted parameter and prove no ordering. Novaro and Wolf [NW71] are described in the citing literature as using an $l(l+1)$ term in a rotor Hamiltonian; we could not obtain the text, and the public record we found shows no ordering proved there. The Laguerre inequality is classical [PS, CV90, Car09], the Gegenbauer facts are in [Sz, DLMF], the stochastic orders, and the sign-change (cut) criterion for the convex order, in [SS07, MS02, KN63], and the shape-invariant structure behind Lemma A in [CKS95] (the partner identity behind Lemma A is already contained there, and Lemma A is not claimed as new). To our knowledge, in the sources we could consult, the convex order between consecutive blocks, its closed form and the first-order theorem do not appear earlier. We could not consult the full texts of Novaro and Wolf (1971), Kitagawara and Barut (1983, 1984), Ostrovsky (1981, 2001, 2003), Allen and Knight (2002) or Wong (1979); they may contain related statements. We also did not read the full text of the Russian original of [DO72], nor Boyd [Boy13] beyond the digits recorded from it (Section 6).

# Setup {#p5.setup}

We use the setting of paper 1 [P1] (https://doi.org/10.5281/zenodo.24522163) and recall only what is needed. The zero-energy coupling problem is
$$-\Delta\psi=Z\,w\,\psi\quad\text{in }D^{1,2}(\mathbb R^3),\qquad w(x)=\frac{2}{|x|(1+|x|/R)^2},$$
and by Theorem A of [P1] it is the circle-invariant conformal Laplacian of the round $S^4$ of radius $\sqrt R$: its eigenvalues are $Z_N=N(N+1)/(2R)$, $N=1,2,\dots$, the eigenspace of group $N$ is spanned by the states $(n,l,m)$ with $n+l=N$, and the Kelvin inversion is the swap of the poles $P_0$ (the nucleus, $r=0$) and $P_\infty$ ($r=\infty$), acting on group $N$ by $(-1)^{N-1}$.

**The angle.** Write
$$c=\cos\theta=\frac{R-r}{R+r},\qquad z=1-c^2=\sin^2\theta ,$$
with $\theta$ measured from $P_0$. The nucleus is $c=+1$, the far pole is $c=-1$, and Kelvin parity is $c\mapsto-c$. This is the orientation of [P1], used throughout; the check scripts of this paper use the opposite sign, which changes nothing for even quantities.

**The blocks.** Put $k=N-1$. A degree-$k$ harmonic of $S^4$ splits into $\mathrm{SO}(4)$ blocks $j=0,\dots,k$ ([P1], Section 3.2, Step 5(i)); block $j$ is spanned by the functions $\sin^j\theta\,C_j(c)\,\eta(\omega)$ with $\eta$ a harmonic of degree $j$ on $S^3$ ($h$ in [P1]). For block $j$ set
$$\lambda_j=j+\tfrac32,\qquad \tau_j=(k-j)(k+j+3),\qquad C_j=C^{(\lambda_j)}_{k-j}(c),$$
that is, the Gegenbauer index, the eigenvalue $\tau_j=d(d+2\lambda_j)$ with $d=k-j$ the degree of $C_j$, and the block polynomial. Fixing $j$ we also write $e=j+2=\lambda_j+\tfrac12$, the exponent of the weight $z^{e}$ of the derivative block.
The volume element is $\sin^3\theta\,d\theta\,d\omega=z\,dc\,d\omega$, so the density of $c$ in block $j$ is the even polynomial of degree $2k+2$
$$p_j(c)=\frac{z^{\,j+1}\,C_j(c)^2}{K_j},\qquad K_j=\int_{-1}^1 z^{\,j+1}C_j^2\,dc. \tag{1}$$
A state is phase-free (invariant under the circle) exactly for $j=2l$ ([P1], Theorem A(3)); the phase-free block $l$ of group $N$ has density $p_{2l}$ (for $l=0,1,\dots$, so $p_{2l}$ is the density of the $l$-th phase-free block), and $c=\cos\theta$ is the natural variable because $\psi_{Nlm}$ depends on $r$ through $c$ only ([P1], Section 3.3(b)). We write $F_l(t)=\Pr_{p_{2l}}(|c|\le t)$, $G_j(t)=\Pr_{p_j}(c\le t)$, and $\mathbb E_p$ for expectation.

**Notation.** Paper 2 is self-contained. It keeps from [P1] the symbols $R$, $w$, $Z$, $Z_N$, $N$, $n$, $l$, $m$, $k$, $j$, $\theta$, $P_0$, $P_\infty$, $C^{(\lambda)}_d$ and $\mathrm{Scal}$. Table 1 lists the symbols of this paper that clash with a symbol of [P1].

| Symbol | Meaning here | Meaning in [P1] |
|---|---|---|
| $z=1-c^2$ | weight variable | $s$ is the root of $s(s+1)=2ZR$ |
| $\lambda_j$, $\tau_j$ | Gegenbauer index and eigenvalue of block $j$ | $\lambda=2Z$; $\mu$ is the spindle parameter |
| $F_l$, $G_j$, $H_l$ | distribution functions and their integral | $F$ is $\Omega\Phi$; $H^1_T$ is a Sobolev space |
| $u_j$, $D_j$ | mean distance, its difference | $U_\mu$ is the Demkov-Ostrovsky potential |
| $\mathbb E_p$ | expectation under the law $p$ | $E_N$ is the eigenspace |
| $\eta$ | harmonic of degree $j$ on $S^3$ | $h_j$ |
| $g$, $h$ | weight change and its Kelvin-even weakening | $g$ is a metric, $h$ a harmonic |
| $\epsilon$ | size of the perturbation | a local cutoff in one proof |
| $c$ | $\cos\theta$; $c_0$ a threshold coupling | $c_\pm$, $c_0$, $c_1$ are constants |

*Table 1. Symbols of this paper that clash with a symbol of [P1].*

![The densities $2p_{2l}(|c|)$ of $|c|$ on $[0,1]$ in the phase-free blocks $l=0,1,2$ of the group $N=5$, from the exact polynomials. The block with small $l$ puts its weight close to the poles ($|c|$ near 1), the block with large $l$ near the equator ($|c|=0$): a larger $l$ moves the weight towards the equator.](figures/fig1_densities.png)

**Three facts on Gegenbauer polynomials and one identity** (DLMF 18.8, 18.9.19; Szegő, *Orthogonal Polynomials*, Thm 3.3.1), used for every $j<k$ with $C=C_j$, $e=j+2$, $\tau=\tau_j$, $K=K_j$:

- (G1) $zC''-2ec\,C'+\tau\,C=0$, equivalently $(z^eC')'=-\tau z^{e-1}C$.
- (G2) $C_j'=2\lambda_jC^{(\lambda_j+1)}_{k-j-1}$, hence $p_{j+1}=z^{e}C_j'^2/(\tau K)$ and $p_j=z^{e-1}C_j^2/K$.
- (G3) $\int z^eC'^2\,dc=\tau K$ (integrate by parts; the boundary term vanishes because $e>0$).
- (G4) (arithmetic, not a Gegenbauer fact) $\tau_j+e(e-1)=(k+1)(k+2)=N(N+1)=2RZ_N$ for every $j$.

The zeros of $C_j$ are real, simple and in $(-1,1)$.

# Why first-order dominance fails {#p5.first-order-fails}

The first-order stochastic order between consecutive blocks fails: their distribution functions cross. For the phase-free pairs, $F_{l+1}-F_l$ changes sign $k-2l-2$ times on $(0,1)$ in the exact count for $N\le12$ (observed, not explained: Lemma A below explains each one-step difference $G_j-G_{j+1}$, not this count). Dominance fails for at least one pair in every group $4\le N\le16$ (exact check, in rational arithmetic on the differences of the distribution functions; the first counterexample is $N=4$, $l=0$), and for the last pair it can hold when $k$ is even (for example $N=5$, $l=1$). The Karlin-Novikoff single-crossing criterion therefore does not apply, and the order that does hold is the convex one.

**Lemma A.** For $0\le j<k$,
$$G_j(t)-G_{j+1}(t)=-\frac{1}{\tau_jK_j}\,z(t)^{j+2}\,C_j(t)\,C_j'(t).$$

The right side is $-z^{e}C\,C'$ up to a constant and changes sign $2(k-j)-1$ times. Proof: Appendix A.1. This is the shape-invariant partner structure of the Pöschl-Teller chain (Cooper, Khare and Sukhatme 1995): the density difference of a state and its partner is an exact derivative.

![Why first-order dominance fails and the convex order holds, for $N=6$ and the pairs $l=0,1$. Left: $F_{l+1}(t)-F_l(t)$ changes sign (the small negative parts are shaded and magnified in the lower panel), so $p_{2l+2}$ does not dominate $p_{2l}$ in the ordinary sense. Right: $H_l(t)=\int_t^1(F_{l+1}-F_l)\,ds$ stays positive on $(0,1)$, which is the convex order and is what the tie-break needs. The negative parts are small: at their deepest, $0.004$ below zero for $l=0$ and $0.017$ for $l=1$.](figures/fig2_dominance_vs_convex.png)

# The convex order {#p5.convex-order}

## Statements {#p5.stmt-convex}

Recall that $X$ precedes $Y$ in the increasing convex order if $\mathbb E\,h(X)\le\mathbb E\,h(Y)$ for every increasing convex $h$, equivalently $\mathbb E(X-t)^+\le\mathbb E(Y-t)^+$ for all $t$ (Shaked and Shanthikumar 2007, Thm 4.A.2); in the convex order the same holds for all convex $h$, equivalently $\mathbb E|X-t|\le\mathbb E|Y-t|$ for all $t$ and equal means (Thm 3.A.1).

In this paper "strict" always means the pointwise statement $D>0$ on the open interval $(-1,1)$ (for Theorem 2, $H_l>0$ on $[0,1)$).

**Theorem 1 (strict convex order of consecutive blocks).** For $k\ge1$ and $0\le j<k$, put $D_j(t)=\mathbb E_{p_j}|c-t|-\mathbb E_{p_{j+1}}|c-t|$. Then
$$D_j(t)=\frac{j+2}{(k+1)(k+2)\,\tau_j\,K_j}\;z(t)^{\,j+3}\,Q_j(t),\qquad Q_j=C_j'^2-\frac{j+1}{j+2}\,C_jC_j'',\tag{4}$$
and $D_j>0$ on $(-1,1)$. In general, for $\lambda>\tfrac12$ and $d\ge1$ the normalised laws $z^{\lambda-1/2}(C^{(\lambda)}_d)^2$ and $z^{\lambda+1/2}(C^{(\lambda+1)}_{d-1})^2$ are in strict convex order, with the analogous closed form (constant $e/((\tau+e(e-1))\tau K)$, $e=\lambda+\tfrac12$, $\tau=d(d+2\lambda)$). The statement is false at $\lambda=\tfrac12$ in the strict pointwise sense.

**Theorem 2 (increasing convex order inside a group).** For every $N\ge3$ and $0\le l<\lfloor(N-1)/2\rfloor$, put $H_l(t)=\int_t^1(F_{l+1}-F_l)\,ds$. Then $H_l=D_{2l}+D_{2l+1}$, that is
$$H_l=z^{2l+3}\!\left[\frac{2l+2}{N(N+1)\,\tau_{2l}K_{2l}}\,Q_{2l}+\frac{(2l+3)\,z}{N(N+1)\,\tau_{2l+1}K_{2l+1}}\,Q_{2l+1}\right],\tag{5}$$
with $\tau_{2l}=(k-2l)(k+2l+3)$, $\tau_{2l+1}=(k-2l-1)(k+2l+4)$, and $H_l>0$ on $[0,1)$, $H_l(1)=0$. Hence $\mathbb E_{p_{2l}}[h(|c|)]>\mathbb E_{p_{2l+2}}[h(|c|)]$ for every bounded increasing convex $h$ on $[0,1]$ that is not constant on $[0,1)$. In particular the mean of $|c|$ strictly decreases in $l$. The number of sign changes of $F_{l+1}-F_l$ is observed, not proved.

Theorem 2 is the convex-order statement of the paper. The physics enters through Theorem 3.

## Proof sketch {#p5.proof-convex}

Each step is proved in the appendix.

1. **Reduction to a mean distance** (A.1). For $0\le t\le1$ and an even law, $\mathbb E(|c|-t)^+=u(t)-t$ with $u(t)=\mathbb E|c-t|$, so $H_l=u_{2l}-u_{2l+2}=D_{2l}+D_{2l+1}$, and $D_j=2\int_{-1}^t(G_j-G_{j+1})$.
2. **Lemma A**: the integrand is $-(\tau_jK_j)^{-1}z^{e}CC'$ (closed form of the CDF difference).
3. **Lemma B** (the identity for the mean distance): with $\Gamma_j=z^{e+1}[C'^2-\frac{e-1}{e}CC'']$ one has $\Gamma_j'=-\frac2e\,N(N+1)\,z^{e}CC'$. So the oscillating integrand of step 2 is an exact derivative, and integrating it gives (4).
4. **Lemma C** (Laguerre): $C'^2-\gamma CC''>0$ for a real-rooted polynomial $C$ with simple zeros and $0<\gamma\le1$; here $\gamma=(j+1)/(j+2)\in[\tfrac12,1)$. This is a weak form of the classical Laguerre inequality (Laguerre; Pólya and Szegő II, Part V; Csordas and Varga 1990; arXiv:0911.1122) and is proved in A.2 with strictness at simple zeros.
5. **Assembly.** (4) and Lemma C give $D_j>0$ on $(-1,1)$; (2) and (3) give $H_l=D_{2l}+D_{2l+1}>0$ on $[0,1)$. For $j=2l+1\le k-1$ both $D_{2l}$ and $D_{2l+1}$ have $k-j\ge1$.
6. **Strictness for increasing convex $h$** (A.3), by Tonelli.

![The function $Q_j=C_j'^2-\frac{j+1}{j+2}\,C_jC_j''$ of (4), divided by its value $Q_j(0)$, on $0\le|c|\le1$ for three blocks $(k,j)=(4,1),(7,2),(10,4)$. It stays strictly positive (log scale), also at the zeros of $C_j$ (circles on the lower edge). Exact positivity for every $N\le30$ is checked by `tiebreak_proof_check.py`.](figures/fig3_Q_positive.png)

# The first-order tie-break {#p5.tie-break}

**Theorem 3.** Fix $R>0$ and a bounded measurable $g$ on $[-1,1]$. For $|\epsilon|\sup|g|<1$ let $Z_{(N,l)}(\epsilon)$ be the eigenvalue of
$$-\Delta\psi=Z\,w\,(1+\epsilon\,g(c))\,\psi\quad\text{in }D^{1,2}(\mathbb R^3),\qquad c=\frac{R-r}{R+r}$$
(nucleus at $c=+1$) in the sector $l$ that equals $Z_N=N(N+1)/(2R)$ at $\epsilon=0$. Then $Z_{(N,l)}$ is analytic near $0$ and
$$Z_{(N,l)}'(0)=Z_N\,\mathbb E_{p_{2l}}[h],\qquad h(x)=-\tfrac12\bigl(g(x)+g(-x)\bigr).$$
If $h$ is increasing and convex on $[0,1]$ and not constant on $[0,1)$, then $Z_{(N,l+1)}'(0)<Z_{(N,l)}'(0)$ for all $N\ge3$ and $0\le l<\lfloor(N-1)/2\rfloor$. Hence for each $N$ there is $\epsilon_N>0$ with lower $n$ first in group $N$ for $0<\epsilon<\epsilon_N$. Equivalently, the Kelvin-even part of $-g$ is convex on $[-1,1]$ and not constant on $(-1,1)$. Only the Kelvin-even part of $g$ acts. $\epsilon_N$ depends on $N$; no uniformity is claimed.

**Proof sketch.** (i) A bounded weight change is form-bounded (Hardy, [P1] (3.1)), so the pencil is analytic in $\epsilon$ (Kato [Kat76]). (ii) $g$ depends on $r$ alone, so $l$ stays a good quantum number, the $m$ degeneracy is exact, and each radial eigenvalue near $Z_N$ is simple: non-degenerate perturbation theory sector by sector. (iii) Hellmann-Feynman for the pencil gives $Z_{(N,l)}'(0)=-Z_N\int wg\psi^2/\int w\psi^2$. (iv) The change of variables $c=(R-r)/(R+r)$ maps the radial density to $p_{2l}$, so the derivative is $-Z_N\mathbb E_{p_{2l}}[g]=Z_N\mathbb E_{p_{2l}}[h]$, the odd part dropping out because $p_{2l}$ is even. (v) Signs: a weaker field needs a larger coupling, so $Z_{(N,l)}$ rises by $\epsilon Z_N\mathbb E_{p_{2l}}[h]$; by Theorem 2 the shift is smaller for larger $l$, so the level with larger $l$ (smaller $n=N-l$) reaches threshold first. Details in A.4.

# The Thomas-Fermi field {#p5.tf-field}

**The nearest Tietz field.** The Thomas-Fermi screening function $\varphi_{\rm TF}(\xi)$ solves $\varphi''=\varphi^{3/2}/\sqrt\xi$ with $\varphi(0)=1$ and $\varphi(\infty)=0$, where $\xi$ is the Thomas-Fermi variable. Its slope at the origin is $B=-1.5880710226113753127\ldots$ (value recorded from Boyd [Boy13], computed to 60 digits by `tf_accurate.py`, agreeing with all 38 digits recorded from [Boy13], and confirmed to 35 digits by the interval sandwich of Appendix A.5). The function was introduced by Fermi [Fer28] and first tabulated accurately by Baker [Bak30]; Tietz [Ti55] approximated it by $(1+a\xi)^{-2}$, which is the weight $w$ of Section 2 with $r/R=a\xi$. We take for $(a,A_{\rm T})$ the pair that minimises, over $a>0$ and $A>0$ jointly, the weighted $L^2$ distance
$$\int_0^\infty\xi\,\bigl(\varphi_{\rm TF}(\xi)-A\,(1+a\xi)^{-2}\bigr)^2\,d\xi$$
(weight $\xi$, the whole half line, the amplitude free; `tietz_extremal.py`, recomputed in `tf_path_check.py`). The minimiser is $a=0.536795$, $A_{\rm T}=1.012783$, and the relative $L^2$ distance of the two functions in this norm is $0.071$. We write $R_{\rm T}$ for the corresponding radius. The radius is fitted, not derived. The amplitude $A_{\rm T}$ multiplies the coupling in every group alike and drops out of every sign, but $a$ does not: $g_{\rm TF}$, $h_{\rm TF}$ and every $S_{(N,l)}$ depend on it, so Corollary 4 is a statement for this one fitted $a$. With $f_{\rm TF}(r)=\varphi_{\rm TF}(r/(aR_{\rm T}))$ the weight change is
$$g_{\rm TF}=\frac{f_{\rm TF}}{A_{\rm T}(1+r/R_{\rm T})^{-2}}-1,\qquad c=\frac{R_{\rm T}-r}{R_{\rm T}+r},$$
and $h_{\rm TF}(x)=-\tfrac12\bigl(g_{\rm TF}(x)+g_{\rm TF}(-x)\bigr)$ is its Kelvin-even weakening. Its value is about $0.5$ near $|c|=1$, where the Thomas-Fermi screening has died away while the Tietz weight has not. It is the far side that drives this: there the Thomas-Fermi tail decays like $\xi^{-3}$ and the Tietz tail like $\xi^{-2}$, which gives about $0.5$, while the nucleus side adds only about $(1-1/A_{\rm T})/2$, roughly $0.006$.

**Normalisation of the couplings.** In Proposition 5 and below, a coupling $Z$ is written in units $R_{\rm T}=1$, that is $\hat Z=ZR_{\rm T}$, so that $\hat Z_N=N(N+1)/2$ at $\epsilon=0$; the hat is dropped. The family $w(1+\epsilon g_{\rm TF})$ divides $A_{\rm T}$ out: at $\epsilon=1$ its weight is $w\,f_{\rm TF}/(A_{\rm T}(1+r/R_{\rm T})^{-2})$. The log values $Z=0.99927$ ($N=1$, $l=0$, $\epsilon=1$) and the threshold coupling $c_0=1.05926488$ of the Thomas-Fermi problem $u''=-c\,\varphi_{\rm TF}(\xi)u/\xi$ in the Thomas-Fermi variable (unit $b=1$) are the same quantity in two normalisations, $c_0=2a\hat Z/A_{\rm T}$ (with $a=0.5367947$, $A_{\rm T}=1.0127834$; `tf_path_check.py` prints both, and the value $c_0$ is the one of `tietz_extremal.py`).

![The Kelvin-even weakening $h_{\rm TF}(|c|)$ of the nearest Tietz field by the Thomas-Fermi field (left), and its second derivative against $u=1-|c|$ on a log axis (right). Left: $h_{\rm TF}$ is about $-0.03$ at the equator and rises to about $0.5$ at the poles. Right: $h_{\rm TF}''$ is positive except in the shaded sliver at the smallest $u$, which Section 6 and Appendix A.6 discuss.](figures/fig4_hTF.png)

**What $h_{\rm TF}$ is not.** $h_{\rm TF}$ is increasing on a fine grid, but it is not convex on all of $[0,1]$. Near $|c|=1$ the Kelvin-even part mixes the nucleus ($r\to0$) with the far tail ($r\to\infty$). The nucleus side carries the cusp term $\tfrac43\xi^{3/2}$ of $\varphi_{\rm TF}$, and it makes $h_{\rm TF}''$ negative on an interval of length about $10^{-7}$ at $|c|=1$ (Appendix A.6). This is real Thomas-Fermi behaviour, not an artefact. It does not matter here, for three reasons. The weight: the share of $H_0$ on $1-t<10^{-3}$ is $4.6\cdot10^{-7}$ at $N=30$, and the probability of $|c|>1-10^{-7}$ under $p_0$, the interval that contains the concave sliver, is $1.4\cdot10^{-10}$ at $N=30$ (two different quantities on two different intervals). The size: $h_{\rm TF}'$ decreases by about $3\cdot10^{-4}$ over the concave interval. The physics: $1-|c|=10^{-7}$ is $r\approx5\cdot10^{-8}R_{\rm T}$ on the nucleus side and $r\approx2\cdot10^{7}R_{\rm T}$ on the far side, where the Thomas-Fermi model describes no atom. The statement "$h_{\rm TF}$ is increasing and convex on $[0,1]$" is therefore not made, and Theorem 3 is not applied to $g_{\rm TF}$ as a whole. The first formula of Theorem 3 needs only that $g$ is bounded, so the first-order slopes are still $Z_N\mathbb E_{p_{2l}}[h_{\rm TF}]$, and Corollary 4 is a direct statement about them.

**Corollary 4 (Thomas-Fermi, certified first order).** Let $S_{(N,l)}=\mathbb E_{p_{2l}}[h_{\rm TF}]-\mathbb E_{p_{2l+2}}[h_{\rm TF}]$. Then $S_{(N,l)}>0$ for every $3\le N\le17$ and every $0\le l<\lfloor(N-1)/2\rfloor$, certified by rigorous interval quadrature of the direct form $S_{(N,l)}=2\int_0^1h_{\rm TF}(p_{2l}-p_{2l+2})\,dx$. The smallest certified lower bound is $0.00685$ ($N=17$, $l=7$). Consequently $Z'_{(N,l+1)}(0)<Z'_{(N,l)}(0)$ for the perturbation $g_{\rm TF}$: in every group with $N\le17$ the first-order shift puts lower $n$ first, and $N\le8$ covers every known element. For $18\le N\le30$ the estimates of $S_{(N,l)}$ are positive for every pair, but certification is not claimed: the enclosures widen with $N$ (Appendix A.5). Monotonicity of $h_{\rm TF}$ is not claimed as certified; it is a numerical observation.

*Proof (computation).* Enclose $\varphi_{\rm TF}$ by a comparison sandwich in interval arithmetic on $|c|\le1-u_0$, $u_0=10^{-3}$, which gives interval enclosures of $h_{\rm TF}$ on 695 pieces. Integrate $2h_{\rm TF}(p_{2l}-p_{2l+2})$ over these pieces by the midpoint rule with its rigorous remainder, the densities being exact rational polynomials. On $|c|>1-u_0$ use only a crude bound for $h_{\rm TF}$ against the exact mass of the densities there. The sum of the lower ends is the certified bound. A second form obtained by integration by parts agrees with it for all 210 pairs with $N\le30$, but its enclosure of $h_{\rm TF}''$ is too wide to carry a proof, so it serves only as a looser cross-check. Details are in Appendix A.5; the code is code/tf_certificate.py (13.5 minutes for $N\le30$, one process).

**Proposition 5 (non-perturbative path; numerical).** For the family of weights $w\,(1+\epsilon\,g_{\rm TF})$, $\epsilon\in[0,1]$ ($\epsilon=0$ the Tietz field, $\epsilon=1$ the Thomas-Fermi field), no level crosses inside a group: $Z_{(N,l+1)}(\epsilon)<Z_{(N,l)}(\epsilon)$ for all $3\le N\le12$, all $0\le l<\lfloor(N-1)/2\rfloor$ and all 67 grid values of $\epsilon\in(0,1]$ tested (68 with $\epsilon=0$) (step $0.025$, refined to $0.0025$ where a relative gap fell below $10^{-3}$). The smallest relative gap over the grid is $3.4\cdot10^{-5}$ ($N=11$, $l=4$, $\epsilon=0.0025$). This minimum sits at the first grid point and is set by the grid step, not by the field: at the Tietz end the gap is linear in $\epsilon$, with slope (relative gap divided by $\epsilon$ at $\epsilon=0.0025$) about $1.4\cdot10^{-2}$ for $N=11$, $l=4$, so a finer grid would give a smaller minimum. The gaps grow along the path. For $\epsilon$ below the grid the order rests on Theorem 3 alone, whose window $\epsilon_N$ is not quantified. At $\epsilon=1$ every group is ordered with higher $l$ first, lower $n$ first, which is the Madelung order within a group (for example $N=3$: $l=1$ at $5.88759$, $l=0$ at $6.49265$; $N=12$: $l=5,4,3,2,1,0$ at $76.357,\dots,91.246$). Method: shooting of the radial Pöschl-Teller equation $y''=(\kappa^2-2ZW)y$, where $\kappa=2l+1$, $\rho=r/R_{\rm T}=e^{2\sigma}$ and $W(\sigma)=4\rho\,(1+\epsilon g_{\rm TF})/(1+\rho)^2$ (it is $\operatorname{sech}^2\sigma$ at $\epsilon=0$). This is the radial equation $u''=[l(l+1)/\rho^2-Z\,w_\epsilon]u$ for $u=\rho R$ after the substitution $u=\rho^{1/2}y(\sigma)$, which turns $\tfrac14\rho^{-3/2}(y_{\sigma\sigma}-y)$ into $u''$ and produces $\kappa^2=1+4l(l+1)$; $Z$ is in the units of the normalisation above. We use the accurate Thomas-Fermi function, fourth-order Numerov, Wronskian matching, $|\sigma|\le12$, double precision. The shooting error is about $10^{-9}$ relative, four orders of magnitude below the gap at the smallest grid value of $\epsilon$. Checks: at $\epsilon=1$, $N=1$, $l=0$ the coupling is $Z=0.99927$, which is the threshold coupling $c_0=2aZ/A_{\rm T}=1.05926488$ of the Thomas-Fermi problem in the normalisation of Section 6; at $\epsilon=0$ the values equal $N(N+1)/2$ to $4\cdot10^{-8}$. This is a grid check in double precision, not a certified enclosure, and between grid points it is a numerical statement. **Overlap of groups at $\epsilon=1$.** Within one group the order is always Madelung's, but the groups themselves begin to overlap: the first overlap is between groups 11 and 12, where $N=11$, $l=0$ lies at $76.974$, above $N=12$, $l=5$ at $76.357$ (the groups $N\le10$ and the pair 10 and 11 do not overlap: $N=10$, $l=0$ at $63.915$, $N=11$, $l=5$ at $64.208$). So in the Thomas-Fermi field the first part of the rule, ordering by $n+l$, is no longer exact between groups 11 and 12, and what is claimed here is the order inside a group only. Real atoms use groups up to $N=8$ only (the 118 known elements fill $n+l\le8$), so the overlap lies beyond the periodic table.

![Non-perturbative check along the path $w\,(1+\epsilon g_{\rm TF})$ from the Tietz field ($\epsilon=0$) to the Thomas-Fermi field ($\epsilon=1$): the levels $Z_{(N,l)}(\epsilon)/Z_N$ for the groups $N=3,\dots,8$ and all $l$. Inside each group the levels spread out in the order of Theorem 3 (lowest $l$ highest coupling) and never cross. The smallest gap over the grid is given in Proposition 5.](figures/fig5_path.png)

**Where the Thomas-Fermi field departs from the round sphere.** The Tietz weight and its nodeless state give a round four-sphere exactly. The following lemma says how to measure the departure of another field. Write the equation as $-\Delta\psi=2Zv\psi$ (so $w=2v$) and let $r=|u|^2$ under the Kustaanheimo-Stiefel map [KS65].

**Lemma D.** The metric $\psi(|u|^2)^2\,\delta$ on $\mathbb R^4$ has scalar curvature $\mathrm{Scal}=48\,Zrv/\psi^2$.

*Proof.* On functions of $r=|u|^2$, $\Delta_4=4r\Delta_3$. In dimension 4, $\mathrm{Scal}(F^2\delta)=-6\,\Delta F/F^3$. With $F=\psi(|u|^2)$ this gives $-24\,r\,\Delta_3\psi/\psi^3=48\,Zrv/\psi^2$. $\square$

For the Tietz weight $v=1/(r(1+r/R)^2)$, the nodeless state $\psi=1/(1+r/R)$ and $Z=1/R$ this is the constant $48/R$ (the round sphere of radius $\sqrt R/2$ in the metric $\psi^2\delta$ on $\mathbb R^4$ used here, which is one quarter of the metric of [P1] Theorem A, whose sphere has radius $\sqrt R$: the two differ by the constant factor $2$ in lengths, not in the conformal class). The constancy is special to this state: for the nodal Tietz states the formula holds but is not constant, so $\mathrm{Scal}=48/R$ is a statement about $N=1$ only. The lemma is checked symbolically in code/h2_scalar_check.py.

Applied to the Thomas-Fermi field at its first threshold coupling ($c_0=1.05926488$ in the Thomas-Fermi normalisation of Section 6, $N=1$, $l=0$, $\psi$ normalised to $\psi(0)=1$), with the round sphere of the nearest Tietz field as reference, the ratio $\mathrm{Scal}_{\rm TF}/\mathrm{Scal}_{\rm round}$ is within $5.3$ percent of $1$ on the nucleus hemisphere $c>0$ (largest deviation $0.0526$, at $c=0.60$). On the far-pole hemisphere the ratio minus $1$ is $+0.008$ at $c=-0.5$, $-0.362$ at $c=-0.9$ and $-0.889$ at $c=-0.99$, and the ratio falls to about $0.003$ at the far pole: it is below $-0.5$ only for $c<-0.9$. The nucleus side also departs, weakly, through the $\xi^{3/2}$ cusp. This is a description on a grid, not a derivation: it shows why the weakening $h_{\rm TF}$ grows toward large $|c|$, since the departure sits at the far pole, which Kelvin parity folds onto the large-$|c|$ end.

![Where the Thomas-Fermi field departs from the round sphere: the ratio $\mathrm{Scal}_{\rm TF}/\mathrm{Scal}_{\rm round}$ of the conformal scalar curvatures against $c$ (the sign convention of [P1]). On the nucleus side ($c>0$) the ratio stays within $5.3$ percent of 1 (largest deviation $0.0526$, at $c$ near $0.6$). It falls away only on the far side, below $0.64$ at $c=-0.9$ and to $0.11$ at $c=-0.99$. The nucleus side also departs weakly through the cusp (Figure 4).](figures/fig6_curvature_ratio.png)

# Discussion {#p5.discussion}

**What is proved, and its limits.** The results are stated in the Abstract and proved in Sections 3 to 6. Theorem 3 is the mechanism: any field near the Tietz field with a Kelvin-even weakening that is increasing and convex in $|\cos\theta|$ puts lower $n$ first. It is first order: its window $\epsilon_N$ depends on $N$ and is not bounded below, and no uniformity in $N$ is claimed. Corollary 4 is certified for $N\le17$; for $18\le N\le30$ the estimates are positive but not certified. Monotonicity of $h_{\rm TF}$ is a numerical observation, and $h_{\rm TF}$ is not convex near the poles, so Theorem 3 does not apply to the Thomas-Fermi field as it stands (Section 6). The reach of first order to the Thomas-Fermi field itself rests on Proposition 5, a grid check in double precision for $N\le12$: it is not certified, and between grid points it is a numerical statement. At the Thomas-Fermi end only the order inside a group is claimed, not the order between groups (see Proposition 5), and real atoms stay in the groups $N\le8$.

**Dependence on the Tietz scale.** Corollary 4 is a statement for the one fitted Tietz scale $a=0.536795$ of Section 6, but that scale is not delicate. A scan over $a$ (`a_scan.py`, $N\le8$ only, whereas Corollary 4 at the fitted $a$ covers $N\le17$) shows that for every consecutive pair with $N\le8$ the difference $S_{(N,l)}$ is positive at every tested $a$ in $[0.40,0.65]$: the grid values are $0.40$, $0.45$, $0.50$, $0.536257$, $0.536795$, $0.60$, $0.65$ (step about $0.05$), with a rigorous interval certificate at $a=0.40$, $a=0.65$ and the fitted $a$, and estimates at the other grid values. Nothing is checked between grid points. By estimate the sign also holds on $0.281<a<0.686$, whose edges are bisection estimates; at $a=0.70$ the pair $(7,2)$ fails. The amplitude $A_{\rm T}$ never enters a sign, and the path at $\epsilon=1$ is the Thomas-Fermi field for every $a$.

**What would close these.** (i) A lower bound for $\epsilon_N$ from the second-order remainder of Theorem 3, which would turn the window into a number. (ii) A validated shooting enclosure along the path, which would turn Proposition 5 into a theorem for the $N$ it covers. (iii) For $N\ge18$, a tighter enclosure of the densities (a Taylor expansion at each piece midpoint instead of the Horner form) or an asymptotic argument in $N$; neither was run. (iv) A tighter enclosure of $h_{\rm TF}'$, which would certify monotonicity.

**Relation to earlier work.** Goudsmit and Richards [GR] interpolate between the Thomas-Fermi limit and the hydrogenic limit. Our base point is the Tietz field, which is exactly solvable and has the degenerate groups, and the interpolation is replaced by a first-order perturbation with a proof; the path from the Tietz field to the Thomas-Fermi field in Proposition 5 is the numerical counterpart of their interpolation. Bjerg, Fournais, Hearnshaw and Solovej [BFHS25] show that in the large-$Z$ many-electron limit the angular momentum distribution follows the Thomas-Fermi one and deviates from Madelung's rule; our statements concern a one-electron model field at low $N$ and say nothing about that limit.

**Open questions.** Second order in $\epsilon$; real self-consistent fields in place of the Thomas-Fermi field; the known exceptions to the rule; and a proof for all $N$ of Corollary 4 or of the order along the path.

# Reproducibility {#p5.repro}


All code is in the repository https://github.com/mik-tf/madelung-tie-break (version v1.0). Each script prints its own verdict line. Run the scripts from the directory `code/`, for example `nice -n 19 env OMP_NUM_THREADS=1 python3 tiebreak_proof_check.py`. From the repository root, `python3 code/run_all.py` runs the quick checks and prints one PASS or FAIL line for each; the flags `--certificate` and `--path` add the two long runs. The exact-arithmetic checks use rational numbers (sympy), so they are exact, not tolerance checks.

| Script | What it checks |
|---|---|
| `tiebreak_proof_check.py` | (G1) to (G3), Lemmas A and B, and the closed forms (4) and (5), exactly, for all $N\le30$ |
| `review_tiebreak_check.py` | independent rebuild of (5) from the states $\psi_{nlm}$ in $\mathbb R^3$ for 90 pairs with $N\le20$; $Q_j$ has no real root; part 4: the first-order slope of Theorem 3 by shooting for $N=3$ to $6$ |
| `tiebreak_dominance.py` | the distribution functions cross: sign-change counts for $N\le12$, failure of dominance for $N\le16$ |
| `tf_accurate.py` | $\varphi_{\rm TF}$ to 60 digits; the slope $B$ and the tail constant $F$; the second differences of $h_{\rm TF}''$ near the poles |
| `tf_certificate.py` | Corollary 4: enclosure, quadrature, cross-check. `nice -n 19 env OMP_NUM_THREADS=1 python3 tf_certificate.py 17 1/1000` prints `VERDICT: PASS`; the run with `NMAX = 30` (log: `code/tf_certificate_N30.log`) prints `VERDICT: FAIL` by design, because the certificate does not close for $N\ge18$ |
| `tf_path_check.py` | Proposition 5: zero-energy couplings along the path; `nice -n 19 env OMP_NUM_THREADS=1 python3 tf_path_check.py` ends with a PASS line (log: `code/tf_path_check.log`) |
| `h2_scalar_check.py` | Lemma D, symbolically, in six steps |
| `a_scan.py` | the dependence on the Tietz scale $a$ (Section 7): estimates on the grid, bisection of the edges, and with `--certify` the interval certificates at $a=0.40$ and $a=0.65$. `nice -n 19 env OMP_NUM_THREADS=1 python3 a_scan.py --certify` (log: `code/a_scan.log`; the two certificates take about 45 s and 50 s) |
| `make_figures.py` | all figures, and the automated layout check: no text box overlaps another text box, an arrow, a marker, a line or the figure edge |

Figures and their data:

| Figure | Data from |
|---|---|
| 1, densities $2p_{2l}$ | exact polynomials as in `tiebreak_proof_check.py`, computed in `make_figures.py` |
| 2, $F_{l+1}-F_l$ and $H_l$ | exact polynomials as in `tiebreak_dominance.py` and `tiebreak_proof_check.py`, computed in `make_figures.py` |
| 3, $Q_j$ | exact polynomials as in `tiebreak_proof_check.py`, computed in `make_figures.py` |
| 4, $h_{\rm TF}$ and the right panel | `tf_certificate.py`, `tf_accurate.py` |
| 5, $Z_{(N,l)}(\epsilon)/Z_N$ | `tf_path_check.py` (via `w3_cache.py`) |
| 6, curvature ratio | computed in `make_figures.py` from the zero-energy state of the Thomas-Fermi weight, with `tf_accurate.py` and the threshold coupling $c_0$ of `tf_path_check.py` |

Point counts, error sizes and run times are in Appendix A.7.

# Use of AI {#p5.ai .unnumbered}

Claude (Anthropic) was used in November 2026, under the author's direction, for exploratory computation, derivations, code, checks, literature search and drafting.

# References {#p5.refs .unnumbered}

[BFHS25] Bjerg, Fournais, Hearnshaw and Solovej, "The angular momentum distribution of electrons in large atoms and its deviation from Madelung's rule", arXiv:2506.14227 (2025).

[Bak30] E. B. Baker, "The application of the Fermi-Thomas statistical model to the calculation of potential distribution in positive ions", Phys. Rev. 36, 630-647 (1930). https://doi.org/10.1103/PhysRev.36.630

[Boy13] J. P. Boyd, "Rational Chebyshev series for the Thomas-Fermi function: Endpoint singularities and spectral methods", J. Comput. Appl. Math. 244, 90-101 (2013). https://doi.org/10.1016/j.cam.2012.11.015 (the digits of $\varphi_{\rm TF}'(0)$ are those recorded from this source in the header of `tf_accurate.py`; the full text was not consulted).

[Car09] D. A. Cardon, "Extended Laguerre inequalities and a criterion for real zeros", arXiv:0911.1122 (2009).

[CKS95] F. Cooper, A. Khare and U. Sukhatme, "Supersymmetry and quantum mechanics", Phys. Rep. 251, 267-385 (1995). https://doi.org/10.1016/0370-1573(94)00080-M

[CM50] C. A. Coulson and N. H. March, "Momenta in atoms using the Thomas-Fermi method", Proc. Phys. Soc. A 63, 367-374 (1950). https://doi.org/10.1088/0370-1298/63/4/306

[CV90] G. Csordas and R. S. Varga, "Necessary and sufficient conditions and the Riemann hypothesis", Adv. Appl. Math. 11, 328-357 (1990). https://doi.org/10.1016/0196-8858(90)90013-O

[DLMF] NIST Digital Library of Mathematical Functions, https://dlmf.nist.gov/, Sections 18.8 and 18.9.

[DO71] Yu. N. Demkov and V. N. Ostrovsky, "Internal symmetry of the Maxwell fish-eye problem and the Fock group for the hydrogen atom", Sov. Phys. JETP 33, 1083 (1971).

[DO72] Yu. N. Demkov and V. N. Ostrovsky, "The n+l filling rule in the periodic system and focusing potentials", Sov. Phys. JETP 35, 66 (1972) (title, authors, volume and first page confirmed in OpenAlex; the Russian original, Zh. Eksp. Teor. Fiz. 62, 125, and the last page 69 could not be confirmed; full text not consulted).

[Fer28] E. Fermi, "Eine statistische Methode zur Bestimmung einiger Eigenschaften des Atoms und ihre Anwendung auf die Theorie des periodischen Systems der Elemente", Z. Phys. 48, 73-79 (1928). https://doi.org/10.1007/BF01351576

[GR] S. A. Goudsmit and P. I. Richards, "The order of electron shells in ionized atoms", Proc. Natl. Acad. Sci. USA 51, 664-671 (1964); Brookhaven report BNL-7786. https://doi.org/10.1073/pnas.51.4.664 (open copy: https://www.osti.gov/servlets/purl/4007960).

[Kat76] T. Kato, *Perturbation Theory for Linear Operators*, 2nd ed., Springer (1976).

[KN63] S. Karlin and A. Novikoff, "Generalized convex inequalities", Pacific J. Math. 13, 1251-1279 (1963). https://doi.org/10.2140/pjm.1963.13.1251

[KS65] P. Kustaanheimo and E. Stiefel, "Perturbation theory of Kepler motion based on spinor regularization", J. reine angew. Math. 218, 204-219 (1965). https://doi.org/10.1515/crll.1965.218.204

[MS02] A. Müller and D. Stoyan, *Comparison Methods for Stochastic Models and Risks*, Wiley (2002).

[NK92] T. Negadi and M. Kibler, "A q-deformed Aufbau Prinzip", J. Phys. A 25, L157-L160 (1992); arXiv:quant-ph/0409177. https://doi.org/10.1088/0305-4470/25/4/002

[NW71] O. Novaro and K. B. Wolf, "A model Hamiltonian for the periodic table", Rev. Mex. Fis. 20, 265-268 (1971) (not found in any open index; unverified; full text not consulted).

[P1] Michaël Perreault, "The Demkov-Ostrovsky Madelung model as the circle-invariant conformal Laplacian of the round 4-sphere", Zenodo (2026), https://doi.org/10.5281/zenodo.24522163.

[PS] G. Pólya and G. Szegő, *Problems and Theorems in Analysis II*, Springer (1976), Part V (the part number is not confirmed against an open copy).

[SS07] M. Shaked and J. G. Shanthikumar, *Stochastic Orders*, Springer (2007). https://doi.org/10.1007/978-0-387-34675-5

[Sz] G. Szegő, *Orthogonal Polynomials*, 4th ed., AMS Colloquium Publications 23 (1975).

[Ti55] T. Tietz, "Über eine Approximation der Fermischen Verteilungsfunktion", Ann. Phys. 450, 186-188 (1955). https://doi.org/10.1002/andp.19554500309

Listed but full text not consulted (see the Introduction): Y. Kitagawara and A. O. Barut, J. Phys. B 16, 3305-3327 (1983), https://doi.org/10.1088/0022-3700/16/18/006, and J. Phys. B 17, 4251-4259 (1984), https://doi.org/10.1088/0022-3700/17/21/013; V. N. Ostrovsky, J. Phys. B 14, 4425 (1981), Found. Chem. 3, 145 (2001) and Ann. N. Y. Acad. Sci. 988, 182 (2003); L. Allen and D. Knight, Int. J. Quantum Chem. 90, 80 (2002); W. Wong, J. Chem. Educ. 56, 714 (1979).

# Appendix A: long proofs {#p5.app-a .unnumbered}

## A.1 Lemmas A and B, and the reduction to a mean distance {#p5.app-a1 .unnumbered}

**Reduction.** For an even law $p$ and $0\le t\le1$ we have $(|c|-t)^+=(c-t)^++(-c-t)^+$ (at most one term is nonzero), the two terms have the same mean by symmetry, and $|c-t|=2(c-t)^+-(c-t)$ with $\mathbb E\,c=0$, so
$$\mathbb E_p(|c|-t)^+=2\,\mathbb E_p(c-t)^+=u_p(t)-t,\qquad u_p(t)=\mathbb E_p|c-t|.$$
By Fubini $\int_t^1(1-F(s))\,ds=\mathbb E(|c|-t)^+$. Hence, for $t\in[0,1]$,
$$H_l(t)=\mathbb E_{p_{2l}}(|c|-t)^+-\mathbb E_{p_{2l+2}}(|c|-t)^+=u_{2l}(t)-u_{2l+2}(t)=D_{2l}(t)+D_{2l+1}(t),\tag{2}$$
with $u_j=u_{p_j}$ (the mean distance of $c$ from $t$ under block $j$) and $D_j=u_j-u_{j+1}$. Since $u_j'(t)=2G_j(t)-1$ and $u_j(-1)=1+\mathbb E_{p_j}c=1$ for every $j$,
$$D_j(t)=2\int_{-1}^t(G_j-G_{j+1})\,du .\tag{3}$$
So it suffices to show $D_j>0$ on $(-1,1)$ for every $j<k$.

**Proof of Lemma A.** Both sides vanish at $t=-1$ because $e>0$. Differentiating the right side with (G1), and then (G2):
$$-\frac1{\tau K}\bigl[(z^eC')'C+z^eC'^2\bigr]=-\frac1{\tau K}\bigl[-\tau z^{e-1}C^2+z^eC'^2\bigr]=p_j-p_{j+1},$$
which is the derivative of the left side. $\square$

**Lemma B.** Let $\Gamma=z^{e+1}\bigl[C'^2-\frac{e-1}{e}CC''\bigr]$ for any solution $C$ of (G1) with parameters $e,\tau$. Then
$$\Gamma'=-\frac{2}{e}\,\bigl(\tau+e(e-1)\bigr)\,z^{e}\,C\,C' .$$
For $C=C_j$ the constant is $-\tfrac2eN(N+1)$ by (G4).

*Proof.* By (G1), $zC''=2ec\,C'-\tau C$, so $\Gamma=z^eT$ with $T=zC'^2-2(e-1)c\,CC'+\frac{e-1}{e}\tau\,C^2$. Then $\Gamma'=z^{e-1}[-2ec\,T+zT']$. Compute $zT'$ term by term, using $z'=-2c$ and replacing every $zC''$ by $2ecC'-\tau C$:

- from $zC'^2$: $z(-2cC'^2)+2C'\,z(zC'')=-2cz\,C'^2+z(4ec\,C'^2-2\tau\,CC')$;
- from $-2(e-1)c\,CC'$: $-2(e-1)z\,(CC'+cC'^2)-2(e-1)c\,C(2ecC'-\tau C)$;
- from $\frac{e-1}{e}\tau C^2$: $2z\frac{e-1}{e}\tau\,CC'$.

Collecting,
$$zT'=2ec\,zC'^2+\Bigl[z\Bigl(-2\tau-2(e-1)+2\tfrac{e-1}{e}\tau\Bigr)-4e(e-1)c^2\Bigr]CC'+2(e-1)\tau\,c\,C^2 .$$
Adding $-2ec\,T=-2ec\,zC'^2+4e(e-1)c^2CC'-2(e-1)\tau c\,C^2$, the $C'^2$ terms cancel, the $c^2CC'$ terms cancel, and the $C^2$ terms cancel, leaving $z[-2\tau(1-\frac{e-1}e)-2(e-1)]CC'=-\frac2e(\tau+e(e-1))\,z\,CC'$. Multiply by $z^{e-1}$. $\square$

**Closed form (4).** $\Gamma_j(-1)=0$ because $e+1>0$, so $\int_{-1}^tz^eCC'=-\frac{e}{2N(N+1)}\Gamma_j(t)$. With Lemma A and (3),
$$D_j=-\frac{2}{\tau K}\int_{-1}^tz^eCC'=\frac{j+2}{N(N+1)\,\tau_jK_j}\,\Gamma_j=\frac{j+2}{N(N+1)\tau_jK_j}\,z^{j+3}Q_j,$$
because $e+1=j+3$ and $(e-1)/e=(j+1)/(j+2)$. The general form for $\lambda>\frac12$ is the same computation with $e=\lambda+\frac12$ and (G4) replaced by $\tau+e(e-1)$. Combining (2) with (4) for $j=2l$ and $j=2l+1$ gives (5). $\square$

At $\lambda=\frac12$ one has $e=1$, the factor $(e-1)/e$ is $0$, and $D$ is a positive multiple of $z^{2}C'^2$, which vanishes at the interior zeros of $C'$: strictness fails pointwise there.

## A.2 Lemma C and the Laguerre identity {#p5.app-a2 .unnumbered}

**Lemma C.** Let $C$ be a real polynomial of degree $d\ge1$ with $d$ simple real zeros $x_1,\dots,x_d$, and $0<\gamma\le1$. Then $C'^2-\gamma CC''>0$ on all of $\mathbb R$.

*Proof.* Write $C=\beta\prod_r(x-x_r)$, $\beta\ne0$. Off the zeros, $C'/C=\sum_i1/(x-x_i)$, so
$$C'^2-CC''=-C^2\bigl(C'/C\bigr)'=C^2\sum_i\frac1{(x-x_i)^2}=\sum_i\Bigl(\frac{C}{x-x_i}\Bigr)^2,$$
and both sides are polynomials, so the identity holds everywhere. The right side is positive at a non-zero $x$, and at a zero $x_r$ the term $i=r$ equals $\bigl(\beta\prod_{i\ne r}(x_r-x_i)\bigr)^2=C'(x_r)^2>0$ because the zero is simple. Finally $C'^2-\gamma CC''=(1-\gamma)C'^2+\gamma(C'^2-CC'')>0$. $\square$

This is a weak form of the classical Laguerre inequality; the sharp polynomial form $(d-1)C'^2-dCC''\ge0$ (Csordas and Varga 1990; arXiv:0911.1122) gives $C'^2-\gamma CC''\ge0$ for $0\le\gamma\le d/(d-1)$, so the factor $\gamma=(j+1)/(j+2)<1$ is not a close case. The real-zero hypothesis is needed ($C=x^2+1$ gives $-2$ at $0$) and holds: $C_j$ is orthogonal for a positive weight, so its $k-j$ zeros are real, simple and in $(-1,1)$ (Szegő, Thm 3.3.1). Hence $Q_j>0$ on $\mathbb R$, and by (4) $D_j>0$ on $(-1,1)$, with $D_j(\pm1)=0$. The case $k-j=1$ ($C$ linear) is covered by Lemma C.

## A.3 Strictness for increasing convex $h$ {#p5.app-a3 .unnumbered}

Let $h$ be bounded, increasing and convex on $[0,1]$. Then on $[0,1)$
$$h(x)=h(0)+h'_+(0)\,x+\int_{(0,1)}(x-t)^+\,d\varpi(t),\qquad \varpi\ge0 \text{ a measure}.$$
By Tonelli,
$$\mathbb E_{p_{2l}}[h(|c|)]-\mathbb E_{p_{2l+2}}[h(|c|)]=h'_+(0)\,H_l(0)+\int_{(0,1)}H_l\,d\varpi .$$
Here $H_l>0$ on $[0,1)$, so the right side is positive unless $h'_+(0)=0$ and $\varpi=0$ on $(0,1)$, that is unless $h$ is constant on $[0,1)$. A jump of $h$ at $x=1$ is invisible since $P(|c|=1)=0$, so "not constant on $[0,1)$" is the exact condition. If $h$ is only $C^1$, integration by parts gives the same with $H_l(1)=0$ killing the boundary term.

## A.4 The first-order step in full {#p5.app-a4 .unnumbered}

*Setting.* The problem is the pencil $-\Delta\psi=Zw(1+\epsilon g)\psi$ in $D^{1,2}(\mathbb R^3)$ with $g$ a bounded measurable function of $r$, written through $c=(R-r)/(R+r)$. The eigenvalue is a Rayleigh quotient $Z=q(\psi)/\int w(1+\epsilon g)\psi^2$, and by the Hardy bound of [P1] (3.1), $\int w\psi^2\le C\,q(\psi)$, a bounded weight change is form-bounded. For $|\epsilon|\sup|g|<1$ the family is analytic of type (B) in $\epsilon$ (Kato, *Perturbation Theory for Linear Operators*, Ch. VII).

*Simple eigenvalues.* A function of $r$ alone commutes with $\mathrm{SO}(3)$. So $l$ remains a good quantum number for every $\epsilon$, the $m$ degeneracy is exact, and in each sector $l$ the radial (Pöschl-Teller) eigenvalue near $Z_N$ is simple. Thus no degenerate perturbation theory is needed, and $Z_{(N,l)}(\epsilon)$, $0\le l\le\lfloor(N-1)/2\rfloor$, is analytic near $0$. (Remark, not needed: multiplication by a function of $\theta$ also commutes with the pole $\mathrm{SO}(4)$; the blocks $j$ are inequivalent $\mathrm{SO}(4)$ types, so it acts on block $j$ by the scalar $\mathbb E_{p_j}[g]$, the same scalar.)

*Derivative.* Hellmann-Feynman for the pencil gives
$$Z_{(N,l)}'(0)=-Z_N\,\frac{\int w\,g\,\psi^2}{\int w\,\psi^2}.$$
Under $c=(R-r)/(R+r)$ the radial density $w\,r^2R_{nl}^2\,dr$ becomes $p_{2l}(c)\,dc$ (equivalently, the unitary lift of [P1] Theorem A). Hence
$$Z_{(N,l)}'(0)=-Z_N\,\mathbb E_{p_{2l}}[g]=Z_N\,\mathbb E_{p_{2l}}[h],\qquad h(x)=-\tfrac12(g(x)+g(-x)),\tag{6}$$
with $Z_N=N(N+1)/(2R)>0$ common to all $l$ in the group. Only the even part of $g$ enters because $p_{2l}$ is even; written as a function of $|c|$, $\mathbb E_{p_{2l}}[h]=\mathbb E_{p_{2l}}[h(|c|)]$.

*Signs.* $h=-g_{\rm even}$ is positive where the field is weakened. A weaker field needs a larger coupling, so $Z_{(N,l)}$ rises by $\epsilon Z_N\mathbb E_{p_{2l}}[h]$. If $\mathbb E_{p_{2l+2}}[h]<\mathbb E_{p_{2l}}[h]$ (Theorem 2 and A.3), then $Z_{(N,l+1)}'(0)<Z_{(N,l)}'(0)$, and by analyticity and the strict inequality of derivatives there is $\epsilon_N>0$ with $Z_{(N,l+1)}(\epsilon)<Z_{(N,l)}(\epsilon)$ for $0<\epsilon<\epsilon_N$: the level with the larger $l$, the smaller $n=N-l$, reaches its zero-energy threshold first, which is lower $n$ first (for example $3p$ before $4s$ at $N=4$). Equivalence with convexity of $-g_{\rm even}$ in $c$ on $[-1,1]$: $h(|c|)=-g_{\rm even}(c)$, and an even function is convex and increasing in $|c|$ on $[0,1]$ iff it is convex on $[-1,1]$ (it is then minimal at $0$). $\square$

## A.5 The Thomas-Fermi function: enclosure and error budget {#p5.app-a5 .unnumbered}

This appendix describes the certificate of Corollary 4 (code/tf_certificate.py). Here $\xi$ is the Thomas-Fermi variable and $x=|c|$.

**Enclosure of $\varphi_{\rm TF}$ without pole asymptotics.** In $t=\sqrt\xi$ the equation becomes a regular first-order system at the nucleus, so the square-root singularity disappears. It is integrated by validated Taylor stepping (order 60, Lagrange remainder on a Picard a-priori box) in interval arithmetic at 480 bits. The comparison sandwich uses $B_{\rm lo}=B^\ast-10^{-35}$, whose solution reaches zero, and $B_{\rm hi}=B^\ast+10^{-35}$, whose solution turns up; the true $\varphi_{\rm TF}$ lies between the two solutions on the whole range $\xi\le1.419\cdot10^5$, which contains $|c|\le1-u_0$ with $u_0=10^{-3}$. An enclosure on all of $[0,\infty)$ is neither needed nor realistic: the unstable mode about the tail $144\,\xi^{-3}$ grows like $\xi^{4.77}$.

**The far-side bound $\varphi_{\rm TF}(\xi)\le144\,\xi^{-3}$ (all $\xi>0$).** The pole part of the certificate needs this bound for every $\xi\ge3.72\cdot10^{3}$ up to infinity, beyond the sandwich range. Put $S(\xi)=144\,\xi^{-3}$. It solves the same equation: $S''=1728\,\xi^{-5}$ and $S^{3/2}/\sqrt\xi=1728\,\xi^{-9/2-1/2}=1728\,\xi^{-5}$. Let $\delta=\varphi_{\rm TF}-S$ on $(0,\infty)$. As $\xi\to0$, $S\to\infty$ and $\varphi_{\rm TF}\to1$, so $\delta\to-\infty$; as $\xi\to\infty$ both functions tend to $0$, so $\delta\to0$. Suppose $\delta>0$ somewhere. Then $\delta$ attains a positive maximum at an interior point, where $\delta''\le0$. But there $\varphi_{\rm TF}>S>0$, so $\delta''=(\varphi_{\rm TF}^{3/2}-S^{3/2})/\sqrt\xi>0$, a contradiction. Hence $\delta\le0$, that is $\varphi_{\rm TF}\le144\,\xi^{-3}$ for all $\xi>0$ (this is the maximum principle for the comparison with the exact Sommerfeld solution).

**The function $h_{\rm TF}$.** Interval enclosures of $h_{\rm TF}$, $h_{\rm TF}'$ and $h_{\rm TF}''$ are computed on 695 geometric pieces of $[0,1-u_0]$, each split in four for the quadrature.

**Direct quadrature (primary).** With $T=1-u_0$ and exact rational densities (sympy, Gegenbauer closed form),
$$S_{(N,l)}=2\int_0^T h_{\rm TF}\,(p_{2l}-p_{2l+2})\,dx+S_{\rm pole},$$
the integral by the midpoint rule on the pieces with remainder $(\mathrm{len}^3/24)\sup|R''|$, $R=h_{\rm TF}(p_{2l}-p_{2l+2})$. The remainder uses a bound for $h_{\rm TF}''$, but it only multiplies $\mathrm{len}^3$, so the wide enclosure of $h_{\rm TF}''$ is harmless. For the pole part let $m$ be the exact mass of $|c|$ above $T$ and $[h_c-w,h_c+w]=[0.50031,0.50704]$ the crude range of $h_{\rm TF}$ there, from $1+B_{\rm lo}\xi\le\varphi\le1$ near the nucleus and $0<\varphi\le\min(1,144\,\xi^{-3})$ far away. The near-nucleus bounds hold for all $\xi\ge0$: $\varphi_{\rm TF}''=\varphi_{\rm TF}^{3/2}/\sqrt\xi>0$, so $\varphi_{\rm TF}$ is convex and lies above its tangent at $0$, $\varphi_{\rm TF}\ge1+B\xi$; since $B_{\rm lo}<B$ we get $1+B_{\rm lo}\xi\le1+B\xi$; and $\varphi_{\rm TF}\le1$ because $\varphi_{\rm TF}'$ increases and $\varphi_{\rm TF}\to0$ forces $\varphi_{\rm TF}'\le0$, so $\varphi_{\rm TF}$ is decreasing from $\varphi_{\rm TF}(0)=1$. The far side is $\xi\ge(2-u_0)/(au_0)\approx3.72\cdot10^{3}$ ($|c|>1-u_0$ gives $\rho\le u_0/(2-u_0)=5.0025\cdot10^{-4}$, and $\xi=1/(a\rho)$). Then $S_{\rm pole}\ge h_c(m_{2l}-m_{2l+2})-w(m_{2l}+m_{2l+2})$. The certified bound is the lower end of the integral plus this.

**Rigour of the arithmetic.** All of it is `mpmath.iv` interval arithmetic with directed rounding, so the claim "certified" holds up to the correctness of that directed rounding. Four small effects, each below $10^{-140}$ and far below every margin ($\ge6.8\cdot10^{-3}$), are closed in the code: the sub-pieces use interval enclosures of their exact end points and midpoints, and adjacent pieces share the same end enclosure, so the true pieces abut exactly and nothing is counted twice; the half-width $w$ of the pole range is rounded up around a point $h_c$, so $[h_c-w,h_c+w]$ contains the range exactly; and the enclosure of the Taylor steps asserts that every requested $t$ lies inside the validated steps (here $t\le61<376.7$) instead of clipping silently.

**Cross-check (integration by parts).** $2\int_0^Th_{\rm TF}(p-p')\,dx=h(T)(F_{2l}-F_{2l+2})(T)-h'(T)H_l(T)+\int_0^TH_lh''$. Its enclosure of $h_{\rm TF}''$ is far too wide (plus or minus $2.8$ near $x=0$ where $h_{\rm TF}''=0.049$), so this form does not close; it is relative width 2 to 10 percent at $N=3$ to $8$ and overlaps the quadrature interval for all 210 pairs. For the same reason the enclosure of $h_{\rm TF}'$ is wider than $h_{\rm TF}'$ itself on 21 of the 695 pieces, all in $x\in[0,0.182]$. This is a limit of the enclosure, not a finding that $h_{\rm TF}$ fails to increase (the numerical $h_{\rm TF}'$ is positive there), and it is why monotonicity is not claimed as certified.

**Result.** Smallest certified lower bound of $S_{(N,l)}$ over $l$ (and the widest quadrature interval):

| $N$ | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|
| bound | 0.1051 | 0.1015 | 0.04924 | 0.05693 | 0.02861 | 0.03537 | 0.01897 | 0.02415 |
| width | $1\cdot10^{-6}$ | 0 | $6\cdot10^{-7}$ | $7\cdot10^{-7}$ | $8\cdot10^{-7}$ | $1\cdot10^{-6}$ | $1.2\cdot10^{-6}$ | $1.5\cdot10^{-6}$ |

| $N$ | 11 | 12 | 13 | 14 | 15 | 16 | 17 |
|---|---|---|---|---|---|---|---|
| bound | 0.01369 | 0.01764 | 0.01046 | 0.01354 | 0.008334 | 0.01079 | 0.00685 |
| width | $2.3\cdot10^{-6}$ | $6\cdot10^{-6}$ | $2.4\cdot10^{-5}$ | $1.2\cdot10^{-4}$ | $5.9\cdot10^{-4}$ | $3.1\cdot10^{-3}$ | $1.6\cdot10^{-2}$ |

Every lower bound is positive for $N\le17$; the smallest is $0.00685$ at $N=17$, $l=7$, and the smallest certified fraction of the midpoint estimate is $0.57$ ($N=17$, $l=0$). For $N=18$ to $30$ the certificate fails: $123$ of the $210$ pairs with $N\le30$ have a non-positive lower bound, although all midpoint estimates are positive. The cause is enclosure width: the intervals grow about five-fold per $N$, because the Horner enclosure of the exact degree $2N$ densities over a piece suffers from alternating-coefficient dependency. A cure, not run, is to expand $p_{2l}-p_{2l+2}$ in a Taylor series of order about 20 at each piece midpoint; projected at about 70 minutes for $N=17$ to $30$.

## A.6 The near-pole expansion of $h_{\rm TF}$ and the defect bound {#p5.app-a6 .unnumbered}

Orientation: in this paper the nucleus is $c=+1$. The scripts use the opposite sign, so pole-side statements below are translated to $c=+1$. Put $u=1-|c|$. The Kelvin-even $h_{\rm TF}(1-u)$ mixes the nucleus side, $r/R_{\rm T}\approx u/2$, with the far side, $r/R_{\rm T}\approx2/u$.

At the nucleus, $\varphi_{\rm TF}=1+B\xi+\tfrac43\xi^{3/2}+O(\xi^2)$ (the Baker series of $\varphi''=\varphi^{3/2}/\sqrt\xi$ [Bak30, CM50]); the $\xi^{3/2}$ term comes from the Thomas-Fermi density diverging like $r^{-3/2}$. At infinity, $\varphi_{\rm TF}=144\,\xi^{-3}(1+F\xi^{-\lambda}+\cdots)$ with $F=-13.27$ and $\lambda=(\sqrt{73}-7)/2=0.772$. Combining the two sides,
$$h_{\rm TF}''(1-u)=36.2\,u^{-0.228}-0.444\,u^{-1/2}+O(1).$$
The cusp term wins as $u\to0$, so $h_{\rm TF}$ is concave for $u$ below about $9\cdot10^{-8}$. The 60-digit values (code/tf_accurate.py, second differences with step $u/50$) are $h_{\rm TF}''=67,139,234,343,385,12.8$ at $u=10^{-2},\dots,10^{-7}$ and $-1716$ at $u=10^{-8}$, against the predicted $-2026$. At $u=10^{-9}$ and $10^{-10}$ the computed values are $+141086$ and $+6874$ against predicted $-9960$ and $-37502$: that is numerics (tail evaluation or difference step), and nothing in the certificate uses that region.

Defect bound: the concave set lies inside $|c|>1-u_0$, where Corollary 4 uses only the crude range of A.5, so the expansions are not part of the certificate; they explain Figure 4. The share of $H_0$ on $1-t<10^{-3}$ is $3.6\cdot10^{-11}$ at $N=4$ and $4.6\cdot10^{-7}$ at $N=30$ ($3.2\cdot10^{-5}$ at $N=120$), the integral of $H_0$ over $[0,1]$ falls like $1/N^2$ ($0.1$ at $N=4$, $2.2\cdot10^{-3}$ at $N=30$), and the decrease of $h_{\rm TF}'$ over the concave set is about $0.9\,\delta^{1/2}\approx3\cdot10^{-4}$ for $\delta\approx10^{-7}$.

## A.7 Details of the checks {#p5.app-a7 .unnumbered}

*Exact checks.* `tiebreak_proof_check.py`: all $k=N-1\le29$. `review_tiebreak_check.py`: 90 pairs; part 4 (first-order slope by shooting): relative error below $10^{-9}$ (recorded $5\cdot10^{-11}$). `tiebreak_dominance.py`: sign-change counts and dominance failures as in Section 3.

*Thomas-Fermi function.* `tf_accurate.py`: 60-digit Taylor integration joined to the Coulson-March series; $\varphi_{\rm TF}'(0)=-1.5880710226113753127186845094\ldots$ agrees with the value recorded from Boyd [Boy13], $-1.5880710226113753127186845094239501095\ldots$, to all digits listed. The certificate does not rest on this value: its sandwich $B^\ast\pm10^{-35}$ is verified by interval arithmetic. (A shorter string, $-1.58807102261137531271861388$, used in an earlier draft, agreed only to 22 digits and is not used.)

*Certificate.* `tf_certificate.py`: sandwich 27 to 36 s; $N\le8$ in 63 s and $N\le14$ in 137 s (one set of timings, rerun on 2026-11-25, one process, `nice 19`, machine load average about 13), $N\le30$ in 808 s (the logged run).

*Path.* `tf_path_check.py`: 41 values of $\epsilon$ (step $0.025$) plus 27 refined points (step $0.0025$), 68 in all including $\epsilon=0$, $N\le12$, every $l$; run time 677 s. Radial equation $y''=(\kappa^2-2ZW)y$, $\kappa=2l+1$, in $\sigma=\tfrac12\ln\xi$ in the code (this differs from the $\sigma$ of Section 6 by the constant $\tfrac12\ln(1/a)$, which only moves the centre of $W$); two-sided shooting, fourth-order Numerov, Wronskian match at $\sigma=0$, $|\sigma|\le12$, double precision. Shooting error: $\max|Z(h=1/256)-Z(h=1/512)|=1.2\cdot10^{-6}$ absolute ($1.3\cdot10^{-8}$ relative), so about $10^{-9}$ relative for the finer step. At $\epsilon=0$, $\max|Z-N(N+1)/2|=4.3\cdot10^{-8}$. Truncation at $|\sigma|=12$ is estimated at $10^{-10}$ and was not varied. The smallest relative gaps at $\epsilon=0.0025$ by group, $N=3$ to $12$: $2.6,2.5,1.2,1.4,0.72,0.88,0.47,0.60,0.34,0.44$ times $10^{-4}$.

*Curvature.* `h2_scalar_check.py` prints PASS after six symbolic steps. The curvature ratio of Section 6 uses the accurate Thomas-Fermi function and the angle $c=(R-r)/(R+r)$.
