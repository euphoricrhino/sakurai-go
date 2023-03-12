package main

import (
	"flag"
	"fmt"
	"math"
	"math/big"
	"os"
	"path/filepath"
)

/*
	example runs

octave --persist `go run *.go --l=9 --x-start=6.0 --x-step=.001 -x-limit=8 --gamma=10.5 --prec=2000`
octave --persist `go run *.go --l=12 --x-start=8.0 --x-step=.001 -x-limit=10 --gamma=12.75 --prec=2000`
octave --persist `go run *.go --l=23 --x-start=18 --x-step=.001 -x-limit=19 --gamma=20.5 --prec=2000`
octave --persist `go run *.go --l=35 --x-start=10000 --x-step=1 -x-limit=100000 --gamma=909.5 --prec=2000`
*/
func main() {
	var (
		l      int
		prec   uint
		gamma  float64
		xstart float64
		xstep  float64
		xlimit float64
	)
	flag.IntVar(&l, "l", 0, "l")
	flag.UintVar(&prec, "prec", 100, "prec")
	flag.Float64Var(&gamma, "gamma", 0.0, "the square root of 2mVR^2/ħ^2")
	flag.Float64Var(&xstart, "x-start", 0.0, "the start of the kR range")
	flag.Float64Var(&xlimit, "x-limit", 0.0, "the limit of the kR range")
	flag.Float64Var(&xstep, "x-step", 0.0, "the step of kR sampling")

	flag.Parse()

	setPrecOnce(uint(prec))

	if l < 0 {
		panic(fmt.Sprintf("invalid --l: %v", l))
	}

	prev, current := genSphericalBessels(l)

	var (
		x     []float64
		delta []float64
		sigma []float64
		peaks []float64
	)
	xval := xstart
	sigmaFactor := math.Pi * 4 * float64(2*l+1)
	for xval < xlimit {
		x = append(x, xval)
		deltaVal := phaseShift(prev, current, xval, gamma)
		delta = append(delta, deltaVal)
		sigmaVal := math.Sin(deltaVal)
		sigmaVal *= sigmaVal / (xval * xval)
		sigmaVal *= sigmaFactor
		sigma = append(sigma, sigmaVal)
		peaks = append(peaks, sigmaFactor/(xval*xval))
		xval += xstep
	}

	filename := filepath.Join(os.TempDir(), "resonance-scattering.m")
	f, err := os.Create(filename)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	fmt.Fprintf(f, "x=%v;\n", x)
	fmt.Fprintf(f, "delta=%v;\n", delta)
	fmt.Fprintf(f, "sigma=%v;\n", sigma)
	fmt.Fprintf(f, "peaks=%v;\n", peaks)
	fmt.Fprintf(f, "subplot(2,1,1); plot(x,delta,';\\delta_l(k)~kR;'); xlabel('kR'); ylabel('\\delta_l(k)'); set(gca,'fontsize',16);\n")
	fmt.Fprintf(f, "subplot(2,1,2); plot(x,sigma,';\\sigma_l(k)~kR;'); xlabel('kR'); ylabel('\\sigma_l(k)'); set(gca,'fontsize',16);\n")
	fmt.Fprintf(f, "hold on;plot(x,peaks,';4\\pi(2l+1)/k^2;');\n")
	fmt.Fprintf(f, "[resVal,resIdx]=max(sigma);xres=x(resIdx);hold on;plot([xres xres],[min(sigma),resVal],sprintf('--k;k_{res}R=%%f;', xres));\n")

	fmt.Println(filename)
}

// Represents a pair of polynomials from which j_l(x) and n_l(x) can be constructed.
// In particular
// j_l(x) = a(1/x)(sin x/x)+b(1/x)(cos x)
// n_l(x) = b(1/x)(sin x)-a(1/x)(cos x / x)
type polyPair struct {
	a []*big.Float
	b []*big.Float
}

// Evaluates j_l(x).
func (pp *polyPair) evalFirst(x float64) float64 {
	va, vb := pp.evalHelper(x)
	// j_l(x) = va*(sin x)/x+vb*(cos x)
	ans := newFromFloat64(math.Sin(x))
	ans.Quo(ans, newFromFloat64(x))
	ans.Mul(va, ans)
	ans.Add(ans, blankFloat().Mul(vb, newFromFloat64(math.Cos(x))))
	ret, _ := ans.Float64()
	return ret
}

// Evaluates n_l(x).
func (pp *polyPair) evalSecond(x float64) float64 {
	va, vb := pp.evalHelper(x)
	// n_l(x) = vb*(sin x) - va*(cos x)/x
	ans := newFromFloat64(math.Cos(x))
	ans.Quo(ans, newFromFloat64(x))
	ans.Mul(va, ans)
	ans.Neg(ans)
	ans.Add(ans, blankFloat().Mul(vb, newFromFloat64(math.Sin(x))))
	ret, _ := ans.Float64()
	return ret
}

func (pp *polyPair) evalHelper(x float64) (*big.Float, *big.Float) {
	bigx := newFromFloat64(1.0 / x)
	eval := newEvaluator(len(pp.a)-1, bigx)
	va := evalPoly(pp.a, eval)
	vb := evalPoly(pp.b, eval)
	return va, vb
}

func evalPoly(c []*big.Float, eval *evaluator) *big.Float {
	maxPower := len(c) - 1
	ans := blankFloat()
	for k := 0; k < len(c); k++ {
		if c[k] != nil {
			ans.Add(ans, blankFloat().Mul(c[k], eval.pow(maxPower-k)))
		}
	}
	return ans
}

// Generates the spherical bessel functions of degree-l and degree-(l+1).
func genSphericalBessels(l int) (*polyPair, *polyPair) {
	f0 := &polyPair{
		a: []*big.Float{newFromFloat64(1)},
		b: []*big.Float{blankFloat()},
	}
	f1 := &polyPair{
		a: []*big.Float{newFromFloat64(1), nil},
		b: []*big.Float{newFromFloat64(-1), nil},
	}
	if l == 0 {
		return f0, f1
	}
	for k := 2; k <= l; k++ {
		step(f0, f1, k)
		f0, f1 = f1, f0
	}
	return f0, f1
}

// p_k=(2k-1)*p_{k-1}-p_{k-2}, according to the recursion
// f_{k+1}=(2k+1)/x f_k - f_{k-1} for f=j_l, or n_l.
// pprev <-> p_{k-2}, prev <-> p_{k-1}, p_k will be overwriting pprev upon output.
func step(pprev, prev *polyPair, k int) {
	factor := newFromInt(2*k - 1)
	va := blankFloat().Mul(prev.a[0], factor)
	vb := blankFloat().Mul(prev.b[0], factor)
	pprev.a = append([]*big.Float{va, nil}, pprev.a...)
	pprev.b = append([]*big.Float{vb, nil}, pprev.b...)
	for l := 2; l <= k; l += 2 {
		va, vb = pprev.a[l], pprev.b[l]
		va.Neg(va)
		vb.Neg(vb)
		if l < len(prev.a) {
			va.Add(va, blankFloat().Mul(prev.a[l], factor))
			vb.Add(vb, blankFloat().Mul(prev.b[l], factor))
		}
	}
}

// Evaluator for x raised to some power, after construction, all power evaluations can be run in logarithmic time.
type evaluator struct {
	x *big.Float
	// Precomputed all powers of form x^(2^n).
	powers []*big.Float
}

func newEvaluator(maxPower int, x *big.Float) *evaluator {
	bits := 0
	n := maxPower + 1
	for n != 0 {
		n = n >> 1
		bits++
	}
	eval := &evaluator{
		x:      x,
		powers: make([]*big.Float, bits+1),
	}
	if bits == 0 {
		return eval
	}
	eval.powers[1] = blankFloat().Set(x)
	for k := 2; k <= bits; k++ {
		eval.powers[k] = blankFloat().Mul(eval.powers[k-1], eval.powers[k-1])
	}
	return eval
}

func (eval *evaluator) pow(n int) *big.Float {
	ans := newFromFloat64(1.0)
	shift := 1
	for n != 0 {
		if n&1 == 1 {
			ans.Mul(ans, eval.powers[shift])
		}
		shift++
		n = n >> 1
	}
	return ans
}

// Compuates the phase shift delta, given x=kR, and gamma=sqrt(2mVR^2/ħ^2).
func phaseShift(prev, current *polyPair, x float64, gamma float64) float64 {
	// let x = kR, then, y = kappa*R, where kappa*R=sqrt(k^2r^2+gamma^2).
	y := math.Sqrt(x*x + gamma*gamma)
	jlx := current.evalFirst(x)
	jly := current.evalFirst(y)
	nlx := current.evalSecond(x)

	l := len(prev.a)
	vx := float64(l+1) / x
	vy := float64(l+1) / y
	// f_l' = f_{l-1}-(l+1)f_l/x
	jlxp := prev.evalFirst(x) - vx*jlx
	jlyp := prev.evalFirst(y) - vy*jly
	nlxp := prev.evalSecond(x) - vx*nlx

	tan := (x*jlxp*jly - y*jlyp*jlx) / (x*nlxp*jly - y*jlyp*nlx)
	// math.arctan() returns [-pi/2,pi/2], so we have to shift the negative values by pi.
	delta := math.Atan(tan)
	if delta < 0 {
		delta += math.Pi
	}
	return delta
}
