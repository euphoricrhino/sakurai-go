package main

import (
	"flag"
	"fmt"
	"math/big"

	"modern.qm/go/util"
)

func main() {
	var l int

	flag.IntVar(&l, "l", 0, "l")

	flag.Parse()

	if l < 0 {
		panic(fmt.Sprintf("invalid --l: %v", l))
	}

	str := "\\begin{aligned} "
	for m := -l; m <= l; m++ {
		str += formula(l, m)
		str += "\\\\"
	}
	str += "\\end{aligned}"
	util.RenderMath(str, "spherical-harmonics.html")
}

// Latex formula for Y_l^m(theta, phi).
func formula(l, m int) string {
	mSign := 1
	if m < 0 {
		m = -m
		mSign = -1
	}
	denomVals := make([]int, 2*m+1)
	denomVals[0] = 4
	for k := 0; k < 2*m; k++ {
		denomVals[k+1] = l + m - k
	}
	denomIn, denomOut := util.SqrtProducts(denomVals)

	numIn, numOut := util.SqrtProducts([]int{2*l + 1})

	pl := constructLegendre(l)
	plm := pl.assoc(m)

	out := util.BlankRat().SetFrac(numOut, denomOut)
	out.Mul(out, plm.overall)

	in := util.BlankRat().SetFrac(numIn, denomIn)

	// Further simplify the denominator of "in" by canceling common factors with numerator of "out".
	g := util.BlankInt().GCD(nil, nil, out.Num(), in.Denom())
	gg := util.BlankRat().SetFrac(util.BlankInt().Mul(g, g), one)
	in.Mul(in, gg)
	out.Mul(out, util.BlankRat().SetFrac(one, g))

	str := fmt.Sprintf("Y_{l=%v}^{m=%v}(\\theta,\\phi)&=", l, m*mSign)
	if mSign >= 0 && m%2 == 1 {
		str += "-"
	}
	str += fmt.Sprintf("\\frac{%v}{%v}", out.Num(), out.Denom())

	inDenomStr := in.Denom().String()
	if in.Denom().Cmp(one) == 0 {
		inDenomStr = ""
	}
	str += fmt.Sprintf("\\sqrt{\\frac{%v}{%v\\pi}}", in.Num(), inDenomStr)
	if m > 1 {
		str += fmt.Sprintf("e^{%vi\\phi}\\sin^{%v}\\theta", m*mSign, m)
	}
	if m == 1 {
		if mSign == -1 {
			str += "e^{-i\\phi}"
		} else {
			str += "e^{i\\phi}"
		}
		str += "\\sin\\theta"
	}
	str += plm.String()
	return str
}

// Legendre polynomial with rational coefficients.
type legendre struct {
	deg   int
	coeff []*big.Rat
}

func constructLegendre(l int) *legendre {
	if l == 0 {
		return &legendre{
			deg:   0,
			coeff: []*big.Rat{big.NewRat(1, 1)},
		}
	}

	if l == 1 {
		return &legendre{
			deg:   1,
			coeff: []*big.Rat{big.NewRat(0, 1), big.NewRat(1, 1)},
		}
	}

	pprev := constructLegendre(0)
	prev := constructLegendre(1)
	for ll := 2; ll <= l; ll++ {
		cur := &legendre{
			deg:   ll,
			coeff: make([]*big.Rat, ll+1),
		}

		for k := 0; k <= ll; k++ {
			cur.coeff[k] = util.BlankRat()
			if k > 0 {
				cur.coeff[k].Mul(big.NewRat(int64(2*ll-1), int64(ll)), prev.coeff[k-1])
			}
			if k <= ll-2 {
				cur.coeff[k].Sub(cur.coeff[k], util.BlankRat().Mul(big.NewRat(int64(ll-1), int64(ll)), pprev.coeff[k]))
			}
		}
		pprev = prev
		prev = cur
	}
	return prev
}

// Associated Legendre function (polynomial in cos\theta) with integer coefficients.
type assocLegendre struct {
	deg     int
	overall *big.Rat
	coeff   []*big.Int
}

// Constructs the associated Legendre function from Legendre polynomial.
func (p *legendre) assoc(m int) *assocLegendre {
	if m < 0 || m > p.deg {
		panic(fmt.Sprintf("unexpected derivative call: deg=%v,m=%v", p.deg, m))
	}
	// Common denominator is gcd(n!,2^n), see http://oeis.org/A060818.
	a := util.Factorial(p.deg)
	b := util.BlankInt().Exp(big.NewInt(2), big.NewInt(int64(p.deg)), nil)
	overallDenom := util.BlankInt().GCD(nil, nil, a, b)

	ret := &assocLegendre{
		deg:   p.deg - m,
		coeff: make([]*big.Int, p.deg-m+1),
	}
	for k := m; k <= p.deg; k++ {
		// Scale by common denominator.
		r := util.BlankRat().Mul(p.coeff[k], util.BlankRat().SetFrac(overallDenom, big.NewInt(1)))
		if r.Denom().Cmp(one) != 0 {
			panic("overall denom not working")
		}
		ret.coeff[k-m] = r.Num()
		for j := 0; j < m; j++ {
			ret.coeff[k-m].Mul(ret.coeff[k-m], big.NewInt(int64(k-j)))
		}
	}
	// Calculate the GCD of all the coefficients.
	g := util.BlankInt().Set(ret.coeff[0])
	for k := 1; k <= ret.deg; k++ {
		g.GCD(nil, nil, g, ret.coeff[k])
	}
	for k := 0; k <= ret.deg; k++ {
		ret.coeff[k].Div(ret.coeff[k], g)
	}
	ret.overall = util.BlankRat().SetFrac(g, overallDenom)
	return ret
}

var (
	one = big.NewInt(1)
)

func (p *assocLegendre) String() string {
	if p.deg == 0 {
		return ""
	}
	str := ""
	first := true
	terms := 0
	for k := 0; k <= p.deg; k++ {
		c := p.coeff[k]
		sgn := c.Sign()
		if sgn == 0 {
			continue
		}
		terms++
		if sgn < 0 {
			str += "-"
			first = false
		} else if !first {
			str += "+"
		} else {
			first = false
		}

		abs := util.BlankInt().Abs(c)
		isOne := abs.Cmp(one) == 0
		if !isOne || k == 0 {
			str += abs.String()
		}
		if k > 1 {
			str += fmt.Sprintf("\\cos^{%v}\\theta", k)
		}
		if k == 1 {
			str += "\\cos\\theta"
		}
	}
	if terms > 1 {
		return fmt.Sprintf("\\left(%v\\right)", str)
	}
	return str
}
