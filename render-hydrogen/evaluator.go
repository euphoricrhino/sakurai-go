package main

import (
	"math/big"

	"github.com/ALTree/bigfloat"
)

// Represents a polynomial.
type polynomial struct {
	coeff []*big.Float
}

// Evaluates the polynomial at x, all powers are computed in logarithmic time.
func (poly *polynomial) eval(x *big.Float) *big.Float {
	pEval := newPowerEvaluator(x, len(poly.coeff)-1)
	ans := blankFloat()
	for k := 0; k < len(poly.coeff); k++ {
		if poly.coeff[k] != nil {
			ans.Add(ans, blankFloat().Mul(poly.coeff[k], pEval.pow(k)))
		}
	}
	return ans
}

// Evaluates wavefunction value for a given point.
type evaluator struct {
	n, l, m int
	// Polynomial part of the radial wavefunction, this is the product of r^l and F(a,c,2r/na_0).
	// radPoly multiplied by e^{-r/na_0} will be the unnormalized radial wavefunction.
	radPoly *polynomial
	// Polynomial part of the angular wavefunction.
	// angularPoly multiplied by sin(theta)^m will be P_l^m(cos(theta)).
	angularPoly *polynomial
}

// Creates new (n,l,m) wavefunction evaluator.
func newEvaluator(cfg *config) *evaluator {
	n, l, m := cfg.N, cfg.L, cfg.M
	eval := &evaluator{
		n:           n,
		l:           l,
		m:           m,
		radPoly:     radialPoly(n, l),
		angularPoly: angular(l, m),
	}
	return eval
}

// Calculates the unnormalized probability density for the (n,l,m) wavefunction for the given point (x,y,z).
func (eval *evaluator) probDensity(x, y, z *big.Float) *big.Float {
	// Radial distance.
	r := blankFloat().Mul(x, x)
	r.Add(r, blankFloat().Mul(y, y))
	r.Add(r, blankFloat().Mul(z, z))
	r.Sqrt(r)

	// cos(theta).
	ct := blankFloat().Quo(z, r)

	rad := eval.radPoly.eval(r)
	// Times e^{-r/na_0}
	exp := bigfloat.Exp(blankFloat().Quo(r, newFromInt(-eval.n)))
	rad.Mul(rad, exp)
	// Square of the radial amplitude.
	rad.Mul(rad, rad)

	// sin(theta)^2m
	sin2 := newFromFloat64(1.0)
	sin2.Sub(sin2, blankFloat().Mul(ct, ct))
	sin2m := newPowerEvaluator(sin2, eval.m).pow(eval.m)

	ang := eval.angularPoly.eval(ct)
	ang.Mul(ang, ang)
	ang.Mul(ang, sin2m)
	return blankFloat().Mul(ang, rad)
}

// Constructs the radial polynomial.
func radialPoly(n, l int) *polynomial {
	a := l + 1 - n
	c := 2*l + 2
	deg := -a
	// radPoly is the product of r^l and F(a,c,2r/na_0).
	radPoly := &polynomial{
		coeff: make([]*big.Float, deg+l+1),
	}
	radPoly.coeff[l] = newFromFloat64(1.0)
	for d := 1; d <= deg; d++ {
		factor := newFromRat((a+d-1)*2, (c+d-1)*(n*d))
		radPoly.coeff[d+l] = blankFloat().Mul(radPoly.coeff[d+l-1], factor)
	}
	return radPoly
}

// Constructs the polynomial part of the angular function.
func angular(l, m int) *polynomial {
	pl := legendre(l)
	poly := &polynomial{
		coeff: make([]*big.Float, len(pl.coeff)-m),
	}
	for k := len(pl.coeff) - 1; k >= m; k -= 2 {
		poly.coeff[k-m] = blankFloat().Set(pl.coeff[k])
		for j := 0; j < m-1; j++ {
			poly.coeff[k-m].Mul(poly.coeff[k-m], newFromInt(k-j))
		}
	}
	return poly
}

// Constructs Legendre polynomial P_l recursively.
func legendre(l int) *polynomial {
	if l == 0 {
		return &polynomial{
			coeff: []*big.Float{newFromFloat64(1.0)},
		}
	}
	if l == 1 {
		return &polynomial{
			coeff: []*big.Float{nil, newFromFloat64(1.0)},
		}
	}

	pprev := legendre(0)
	prev := legendre(1)
	for ll := 2; ll <= l; ll++ {
		cur := &polynomial{
			coeff: make([]*big.Float, ll+1),
		}

		factor1 := newFromRat(2*ll-1, ll)
		factor2 := newFromRat(ll-1, ll)
		// P_l only has powers with the same parity as l.
		for k := ll; k >= 0; k -= 2 {
			cur.coeff[k] = blankFloat()
			if k > 0 {
				cur.coeff[k].Mul(prev.coeff[k-1], factor1)
			}
			if k <= ll-2 {
				cur.coeff[k].Sub(cur.coeff[k], blankFloat().Mul(pprev.coeff[k], factor2))
			}
		}
		pprev = prev
		prev = cur
	}
	return prev
}

// Evaluator for x raised to some power, after construction, all power evaluations can be run in logarithmic time.
type powerEvaluator struct {
	x *big.Float
	// Precomputed all powers of form x^(2^n).
	powers []*big.Float
}

// Constructs power evaluator given the maximum possible power to be computed subsequently.
func newPowerEvaluator(x *big.Float, maxPower int) *powerEvaluator {
	bits := 0
	n := maxPower
	for n != 0 {
		n = n >> 1
		bits++
	}
	pEval := &powerEvaluator{
		x:      x,
		powers: make([]*big.Float, bits+1),
	}
	// powers[0] was never explicitly used by pow() below, so let it remain nil.
	if bits == 0 {
		return pEval
	}
	pEval.powers[1] = blankFloat().Set(x)
	for k := 2; k <= bits; k++ {
		pEval.powers[k] = blankFloat().Mul(pEval.powers[k-1], pEval.powers[k-1])
	}

	return pEval
}

func (pEval *powerEvaluator) pow(n int) *big.Float {
	ans := newFromFloat64(1.0)
	shift := 1
	for n != 0 {
		if n&1 == 1 {
			ans.Mul(ans, pEval.powers[shift])
		}
		shift++
		n = n >> 1
	}
	return ans
}
