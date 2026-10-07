# Day 228 PROVE — two-row specialisation of Thm 2.5; the separator, gauge-normalised; Prop M from Cor pieri

Date 2026-10-07. Scripts: `scripts/day228/tworow.py` (Item 5), `scripts/day228/sep.py` (Item 1).

---------------------------------------------------------------------------------------------------
## A. Two-row specialisation of the two-point formula (Thm 2.5 = FPSAC `thm:2pt`)

**Notation.** $p_\mu=\sum_\lambda X^\lambda_\mu(t)P_\lambda(x;t)$ (Macdonald III.7 convention). $\lambda=(\lambda_1,\lambda_2)$, $\lambda_2\ge1$,
$n=|\lambda|$, $\rho=\lambda-(1,1)$, $m=\lambda_1-\lambda_2$. Class $(x,y)$, $x+y=n$, $x,y\ge1$; $m_{xy}=2$ if $x=y$, else $1$.

**Theorem (two rows, two parts).** For $y\le x$:
$$X^{(\lambda_1,\lambda_2)}_{(x,y)}(t)=\begin{cases}(t-1)\,t^{\lambda_2-1-y}(1+t^{y}) & y<\lambda_2,\\ m_{xy}-(1-t)\,t^{\lambda_2-1} & y=\lambda_2,\\ (t-1)\,t^{\lambda_2-1} & y>\lambda_2.\end{cases}$$
(The class is unordered, and the formula from Thm 2.5 is the same whether one inserts $y$ or $x$; this was checked, see C.)

### Proof (substitution into Thm 2.5)
Inputs: Thm 2.5 (proved, cold recheck Day 226) and the Green dictionary
$\Phi_a(P_\rho;x,y)=(1-t^x)(1-t^y)X^\lambda_{(x,y)}/b_\lambda$, $\lambda=\rho+1^a$, which is in the FPSAC remark and was checked symbolically for $|\lambda|\le8$ on Day 226.

1. **$a=2$.** Only $A=B=1$ occurs. $G_1(w)=P_\rho(1,w;t)$.
   In two variables $P_{\rho}=(x_1x_2)^{\rho_2}P_{(m)}$, and from the symmetrisation definition (with $v_{(m,0)}=1$ for $m\ge1$),
   $P_{(m)}(x_1,x_2)=\frac{x_1^m(x_1-tx_2)-x_2^m(x_2-tx_1)}{x_1-x_2}$. Hence
   $$\hat G(w):=G_1(w)=w^{\rho_2}\frac{(1-tw)-w^m(w-t)}{1-w}\quad(m\ge1),\qquad \hat G(w)=w^{\rho_2}\quad(m=0).$$
   This gives $\hat G(t)=t^{\rho_2}(1+t)$ for $m\ge1$, and $t^{\rho_2}$ for $m=0$.
   (PROVE.md's "$1+w^m+(1-t)\sum$" is correct for $m\ge1$ only. At $m=0$ it would give $2w^{\rho_2}$, which is wrong.)
2. **Thm 2.5 at $a=2$.** $(-1)^a=1$, $t^{-AB}=t^{-1}$. Using $\sum_{j\ge1}(t^{-j}-t^j)w^j=\frac{(1-t^2)\,w}{(t-w)(1-tw)}$, the bracket is
   $$\frac{X}{b_\lambda}=\frac1t[w^{y-1}]\hat G(w)\Bigl(\frac1{(1-t)^2}+\frac{w}{(t-w)(1-tw)}\Bigr)-\frac{\hat G(t)(1+t^{-y})}{1-t^2}.$$
3. **Case $m\ge1$.** Here $\hat G\cdot\frac{w}{(t-w)(1-tw)}=w^{\rho_2+1}\Bigl[\frac1{(1-w)(t-w)}+\frac{w^m}{(1-w)(1-tw)}\Bigr]$,
   because $(1-tw)-w^m(w-t)$ splits against $(t-w)(1-tw)$. Put $u=y-\lambda_2$ and $v=y-\lambda_1$. Then
   $[w^{N}]\frac1{(1-w)(t-w)}=t^{-N-1}[N+1]$ and $[w^N]\frac1{(1-w)(1-tw)}=[N+1]$, so the coefficient is
   $$\tfrac{1}{(1-t)^2}\bigl(\delta_{u0}+(1-t)1_{u\ge1}+t\delta_{v0}-(1-t)1_{v\ge1}\bigr)+t^{-u}[u]1_{u\ge0}+[v]1_{v\ge0}.$$
   Take $y\le x$. Then $y\le n/2<\lambda_1$, so $v<0$ and the $v$-terms vanish. With $b_\lambda=(1-t)^2$:
   - $u<0$: $X=-(1-t)t^{\lambda_2-1}(1+t^{-y})$.
   - $u=0$: $X=\frac1t-(1-t)(t^{\lambda_2-1}+t^{-1})=1-t^{\lambda_2-1}+t^{\lambda_2}$.
   - $u\ge1$: $(1-t)X/b=\frac1t+t^{-u-1}-\frac1t-t^{\lambda_2-1}-t^{\lambda_2-1-y}$. Since $t^{-u-1}=t^{\lambda_2-1-y}$, the two cancel, leaving $X=-(1-t)t^{\lambda_2-1}$.
4. **Case $m=0$** ($\lambda=(L,L)$, $b_\lambda=(1-t)(1-t^2)$, $u=y-L\le0$ for $y\le x$).
   The coefficient is $\delta_{u0}/(1-t)^2$, since $[w^u]\frac{w}{(t-w)(1-tw)}=0$ for $u\le0$.
   - $u<0$: $X=-(1-t)t^{L-1}(1+t^{-y})$, the same as in case 3.
   - $u=0$: $X=\frac{1+t}{t}-(1-t)(t^{L-1}+t^{-1})=2-t^{L-1}+t^L$. ∎

**Sanity, $t=1$.** $X^\lambda_\mu(1)=[m_\lambda]p_\mu$. For $y\le x$ this is $m_{xy}\cdot1_{y=\lambda_2}$ ✓.

**Grade:** proved, by substitution into the proved Thm 2.5 plus the proved dictionary. Computed: 155/155 ordered $(\lambda,(x,y))$ with $|\lambda|\le10$
against `green.py` (Kostka–Foulkes × Murnaghan–Nakayama), for both the raw Thm 2.5 evaluation and the closed form.

**Novelty: none claimed.** $\ell(\lambda)=2$ is the range of Morris's formula (LNM 579, 1977, scope per Jing–Liu 2104.04411 p.12; not seen
first-hand). This is the $\ell(\lambda)=2\cap\ell(\mu)=2$ overlap. Its use is as a cross-check for Clio's independent two-row × two-part
derivation, and as a worked instance of Thm 2.5.

---------------------------------------------------------------------------------------------------
## B. The separator, in one normalisation (FPSAC remark "Prior art — graph half")

Gauge freedom on full-merge leads: $L(\lambda)\mapsto L(\lambda)\,a^{\ell(\lambda)}g(|\lambda|)\prod_if(\lambda_i)$. The three factors come from
rescaling generators $e_k\mapsto f(k)e_k$ (which also contributes $1/f(n)$ on the output, absorbed into $g$), from degree-wise
normalisation of a pairing, and from rescaling the expansion parameter.

- The old $I=L(111)/(L(21)L(11))$ is invariant under $f$ and $a$ but **not under $g$**: it changes by $1/g(2)$. So it is **not a separator**.
  Example: HT's $\varphi$-reading $\langle\kappa_{HT}(\lambda),e_n\rangle$ differs from $\star$ by $g(n)=t^{\binom n2}$ modulo $f,a$. That is
  exactly why the "different" value $(t+2)/(t(t+1))$ appeared.
- $J:=L(1^4)L(22)/L(211)^2$ is invariant under all three. Exponent check: parts $m_1$: $4-2\cdot2=0$; $m_2$: $2-2=0$;
  degree: $1+1-2=0$; $\ell$: $4+2-2\cdot3=0$.
- Results (`sep.py`):
  - $\star$ (Thm G): $J=\frac{(t^2+1)(t^3+3t^2+6t+6)}{(t+1)(t^2+t+2)^2}$, with $J(0)=3/2$ and $J(1)=1$.
  - Haglund–Tewari $\langle\kappa,e_n\rangle$: **identical**.
  - Dołęga $\oplus$ column cumulant (lowest $(q-1)$-coefficient of $[g_n]\widetilde H_{\lambda'}$): $J\equiv3/2$.
- **Conclusion.** HT's single-row cumulant is the same graph object (already credited). Dołęga's $\oplus$ is genuinely different for $t\ne0$
  and agrees at $t=0$, consistent with 217e Thm B = DFK Cor 5.18.
- Observation (computed, $n\le4$): Dołęga's $\oplus$-lead is a pure gauge of its $t=0$ value $(\ell-1)!$:
  $L_D=(\ell-1)!(t-1)^{1-\ell}g(n)$ with $g(2)=1$, $g(3)=1/(1+t)$, $g(4)=1/((1+t)(1+t+t^2))$.

Grade: computed (symbolic, exact). The invariance of $J$ is proved by the exponent count above.

---------------------------------------------------------------------------------------------------
## C. Prop M from Cor pieri (FPSAC `prop:M`)

$e_k\star e_r=\sum_y s^yC_{k-y,r-y}e_{k+r-y}e_y$, with $C_{a,b}=(-1)^b\frac{[a+b]}{[a]}\pi_a(b)$ and $\pi_a(b)=[u^b]\prod_{m<a}\frac{1+sut^m}{1+ut^m}$.
At $s=1$, $\pi_a(b)=\delta_{b0}$. The $y=r$ term contributes $r\,e_ke_r$. For $y=r-j$:
$\partial_s\pi_a(j)|_{s=1}=[u^j]\sum_{m<a}\frac{ut^m}{1+ut^m}=(-1)^{j-1}[a]_{t^j}$. Hence
$\partial_sC_{a,j}|_{s=1}=-\frac{[a+j]}{[a]}[a]_{t^j}=\frac{(1-t^{a+j})(t^{aj}-1)}{(1-t^a)(1-t^j)}=L(a,j)$ with $a=k-r+j$. ∎
(sympy check, $a\le4$, $j\le5$: True.) No circularity: the proof of Cor pieri uses Thm DS, the Column Lemma and Thm lin, but not Prop M.
