package main

import (
	"flag"
	"fmt"
	"math/big"

	"modern.qm/go/util"
)

func main() {
	var n int

	flag.IntVar(&n, "n", 0, "n")

	flag.Parse()

	if n < 0 {
		panic(fmt.Sprintf("invalid --n: %v", n))
	}

	str := "\\begin{aligned} "
	str += formulaj(n) + "\\\\"
	str += formulan(n)
	str += "\\end{aligned}"
	util.RenderMath(str, "spherical-bessel.html")
}

// Latex formula for j_l(x).
func formulaj(n int) string {
	if n == 0 {
		return "j_0 &= \\frac{\\sin x}{x}"
	} else if n == 1 {
		return "j_1 &= \\frac{\\sin x}{x^2}-\\frac{\\cos x}{x}"
	}
	a, b := genCoefficients(
		[]*big.Int{big.NewInt(1)},
		[]*big.Int{util.BlankInt()},
		[]*big.Int{big.NewInt(1), nil},
		[]*big.Int{big.NewInt(-1), nil},
		n,
	)
	return fmt.Sprintf("j_{%v}&=", n) + render(a, "\\frac{\\sin x}{x}") + "+" + render(b, "\\cos x")
}

// Latex formula for n_l(x).
func formulan(n int) string {
	if n == 0 {
		return "n_0 &= -\\frac{\\cos x}{x}"
	} else if n == 1 {
		return "n_1 &= -\\frac{\\cos x}{x^2}-\\frac{\\sin x}{x}"
	}
	a, b := genCoefficients(
		[]*big.Int{util.BlankInt()},
		[]*big.Int{big.NewInt(-1)},
		[]*big.Int{big.NewInt(-1), nil},
		[]*big.Int{big.NewInt(-1), nil},
		n,
	)
	return fmt.Sprintf("n_{%v}&=", n) + render(a, "\\sin x") + "+" + render(b, "\\frac{\\cos x}{x}")
}

func genCoefficients(a0, b0, a1, b1 []*big.Int, n int) ([]*big.Int, []*big.Int) {
	for k := 2; k <= n; k++ {
		step(&a0, a1, k)
		step(&b0, b1, k)
		a0, a1 = a1, a0
		b0, b1 = b1, b0
	}
	return a1, b1
}

// p_k=(2k-1)*p_{k-1}-p_{k-2}, according to the recursion
// f_{k+1}=(2k+1)/x f_k - f_{k-1} for f=j_l, or n_l.
// pprev <-> p_{k-2}, prev <-> p_{k-1}, p_k will be overwriting pprev upon output.
func step(pprev *[]*big.Int, prev []*big.Int, k int) {
	factor := big.NewInt(int64(2*k - 1))
	v := util.BlankInt().Mul(prev[0], factor)
	*pprev = append([]*big.Int{v, nil}, *pprev...)
	for l := 2; l <= k; l += 2 {
		v = (*pprev)[l]
		v.Neg(v)
		if l < len(prev) {
			v.Add(v, util.BlankInt().Mul(prev[l], factor))
		}
	}
}

func render(a []*big.Int, term string) string {
	n := len(a) - 1
	str := "\\left("
	for l := 0; l < len(a); l += 2 {
		pow := ""
		if n-l > 1 {
			pow = fmt.Sprintf("^{%v}", n-l)
		}
		sign := ""
		if a[l].Sign() == 0 {
			continue
		}
		if a[l].Sign() < 0 {
			sign = "-"
			a[l].Abs(a[l])
		} else if l > 0 {
			sign = "+"
		}
		if l == len(a)-1 {
			str += fmt.Sprintf("%v%v", sign, a[l])
		} else {
			str += fmt.Sprintf("%v\\frac{%v}{x%v}", sign, a[l], pow)
		}
	}
	str += "\\right)" + term
	return str
}
