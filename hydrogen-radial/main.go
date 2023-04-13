package main

import (
	"flag"
	"fmt"
	"math/big"

	"github.com/euphoricrhino/sakurai-go/util"
)

func main() {
	var n int

	flag.IntVar(&n, "n", 0, "n")

	flag.Parse()

	if n <= 0 {
		panic(fmt.Sprintf("invalid --n: %v", n))
	}

	str := "\\begin{aligned} "
	for l := 0; l < n; l++ {
		str += formula(n, l)
		str += "\\\\"
	}
	str += "\\end{aligned}"
	util.RenderMath(str, "hydrogen-radial.html")
}

var (
	one = big.NewInt(1)
)

// Latex formula for R_{ln}(r).
func formula(n, l int) string {
	p := constructPoly(n, l)
	vals := make([]int, 0, 2*l+1)
	for k := -l; k <= l; k++ {
		vals = append(vals, n+k)
	}
	in, out := util.SqrtProducts(vals)

	// 2^(l+1).
	exp := util.BlankInt().Exp(big.NewInt(2), big.NewInt(int64(l+1)), nil)
	outNum := out.Mul(out, exp)
	outDenom := util.Factorial(2*l + 1)
	outDenom.Mul(outDenom, util.BlankInt().Exp(big.NewInt(int64(n)), big.NewInt(int64(l+2)), nil))
	r := util.BlankRat().SetFrac(outNum, outDenom)

	lTerm := ""
	if l >= 1 {
		lTerm = "\\left(\\frac{Zr}{a_0}\\right)"
	}
	if l > 1 {
		lTerm += fmt.Sprintf("^{%v}", l)
	}

	rootTerm := ""
	if in.Cmp(one) > 0 {
		rootTerm = fmt.Sprintf("\\sqrt{%v}", in)
	}
	overallUp := ""
	if r.Num().Cmp(one) == 0 {
		if rootTerm == "" {
			overallUp = "1"
		} else {
			overallUp = rootTerm
		}
	} else {
		overallUp = r.Num().String() + rootTerm
	}
	overall := ""
	if r.Denom().Cmp(one) == 0 {
		overall = overallUp
	} else {
		overall = fmt.Sprintf("\\frac{%v}{%v}", overallUp, r.Denom())
	}
	expTerm := ""
	if n == 1 {
		expTerm = "e^{-Zr/a_0}"
	} else {
		expTerm = fmt.Sprintf("e^{-Zr/%va_0}", n)
	}
	str := fmt.Sprintf(
		"R_{n=%v,l=%v}(r)&=%v\\left(\\frac{Z}{a_0}\\right)^{3/2}%v%v%v",
		n, l,
		overall,
		lTerm,
		expTerm,
		p,
	)
	return str
}

// Polynomial F(a,c,x), see Sakurai pp203, eq 3.310.
type poly struct {
	deg   int
	coeff []*big.Rat
}

func constructPoly(n, l int) *poly {
	// In language of F(a,c,x), where x=Zr/a0.
	a := l + 1 - n
	c := 2*l + 2
	p := &poly{
		deg:   -a,
		coeff: make([]*big.Rat, -a+1),
	}

	p.coeff[0] = big.NewRat(1, 1)
	for d := 1; d <= p.deg; d++ {
		p.coeff[d] = util.BlankRat().Mul(p.coeff[d-1], big.NewRat(int64(a+d-1)*2, int64(c+d-1)*int64(n*d)))
	}

	return p
}

func (p *poly) String() string {
	if p.deg == 0 {
		return ""
	}
	str := "1"
	const xTerm = "\\left(\\frac{Zr}{a_0}\\right)"
	for d := 1; d <= p.deg; d++ {
		num := p.coeff[d].Num()
		if num.Sign() < 0 {
			str += "-"
		} else {
			str += "+"
		}
		str += fmt.Sprintf("\\frac{%v}{%v}", util.BlankInt().Abs(num), p.coeff[d].Denom())
		if d == 1 {
			str += xTerm
		} else {
			str += fmt.Sprintf("%v^{%v}", xTerm, d)
		}
	}
	return fmt.Sprintf("\\left[%v\\right]", str)
}
