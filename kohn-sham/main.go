package main

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sync"

	"gonum.org/v1/gonum/integrate"
	"gonum.org/v1/gonum/mat"
)

const (
	limit   = 200.0
	samples = 50000
	step    = limit / float64(samples)
)

const (
	zEff  = 2.0 - 5.0/16.0
	zEff2 = zEff * zEff
	zEff3 = zEff2 * zEff
	zEff4 = zEff2 * zEff2
	// One hatree in ev
	hatree = 27.2114
)

// All the immutable radial points in the integral.
var rsamples []float64

func init() {
	for i := 0; i <= samples; i++ {
		rsamples = append(rsamples, step*float64(i))
	}
}

// Generates the data slice for the radial integral, every data point is r^2f(r).
func radialSample(f func(int) float64) []float64 {
	const workers = 200
	var wg sync.WaitGroup
	wg.Add(workers)

	ret := make([]float64, len(rsamples))
	for w := 0; w < workers; w++ {
		go func(w int) {
			for i := w; i < len(rsamples); i += workers {
				r := rsamples[i]
				ret[i] = f(i) * r * r
			}
			wg.Done()
		}(w)
	}
	wg.Wait()
	return ret
}

func radialIntegrate(f func(int) float64) float64 {
	return integrate.Simpsons(rsamples, radialSample(f))
}

func vExt(i int) float64 {
	// r=0 case, ok to return zero, since this will eventually be multiplied by r^2 in the radial integral.
	if i == 0 {
		return 0.0
	}
	return -2.0 / rsamples[i]
}

func uEESlice(n func(int) float64) []float64 {
	ret := make([]float64, len(rsamples))
	for i := range ret {
		integrand := func(ii int) float64 {
			if i == 0 && ii == 0 {
				// Both r and r' are zero, we are divergent here, but safe to return zero since eventually it will integrate with r^2 multiplied in.
				return 0.0
			}
			rg := rsamples[i]
			if rg < rsamples[ii] {
				rg = rsamples[ii]
			}
			return 4.0 * math.Pi * n(ii) / rg
		}
		ret[i] = radialIntegrate(integrand)
	}
	return ret
}

// Functional derivative of U_xc over n.
func fnDerUXC(i int, n func(int) float64) float64 {
	nVal := n(i)
	return epsilonXC(nVal) + derEpsilonXC(nVal)
}

func vKS(i int, uEE, n func(int) float64) float64 {
	return vExt(i) + uEE(i) + fnDerUXC(i, n)
}

func epsilonXC(nVal float64) float64 {
	return epsilonX(nVal) + epsilonC(nVal)
}

const (
	epsilonXA = -1.174
	epsilonXB = 3.0 * math.Pi * math.Pi
)

func epsilonX(nval float64) float64 {
	return epsilonXA * 3.0 / (4.0 * math.Pi) * math.Pow(epsilonXB*nval, 1.0/3.0)
}

const (
	epsilonCA = -0.0233504
	epsilonCB = 0.1018
	epsilonCC = 0.102582
)

func epsilonC(nVal float64) float64 {
	y := math.Pow(4.0*math.Pi*nVal/3.0, 1.0/6.0)
	y2 := y * y
	return epsilonCA * y2 / (y2 + epsilonCB*y + epsilonCC)
}

func derEpsilonXC(nVal float64) float64 {
	return derEpsilonX(nVal) + derEpsilonC(nVal)
}

func derEpsilonX(nVal float64) float64 {
	return epsilonXA * 3.0 * math.Pi / 4.0 * math.Pow(epsilonXB, -2.0/3.0) * math.Pow(nVal, 1.0/3.0)
}

func derEpsilonC(nVal float64) float64 {
	y := math.Pow(4.0*math.Pi*nVal/3.0, 1.0/6.0)
	y2 := y * y
	y3 := y2 * y
	d := y2 + epsilonCB*y + epsilonCC
	return epsilonCA * (epsilonCB/6.0*y3 + epsilonCC/3.0*y2) / (d * d)
}

// Radial functions (leaving out the zEff^3/2 factor, which will be compensated overall at last after the radial integration).
func r1s(i int) float64 {
	return 2.0 * math.Exp(-zEff*rsamples[i])
}

func r2s(i int) float64 {
	r := rsamples[i]
	return (1.0 - zEff*0.5*r) * math.Exp(-0.5*zEff*r) / math.Sqrt(2.0)
}

func r3s(i int) float64 {
	r := rsamples[i]
	return (1.0 - 2.0*zEff*r/3.0 + 2.0*zEff2*r*r/27.0) * math.Exp(-zEff*r/3.0) * 2.0 * math.Sqrt(3) / 9.0
}

// -∇^2[R_1s(r)]/2
func kinR1s(i int) float64 {
	if i == 0 {
		return 0.0
	}
	r := rsamples[i]
	return -(zEff2 - 2.0*zEff/r) * math.Exp(-zEff*r)
}

// -∇^2[R_2s(r)]/2
func kinR2s(i int) float64 {
	if i == 0 {
		return 0.0
	}
	r := rsamples[i]
	return -0.25 * math.Sqrt(2) * (-zEff3*r/8.0 + 1.25*zEff2 - 2.0*zEff/r) * math.Exp(-0.5*zEff*r)
}

// -∇^2[R_1s(r)]/2
func kinR3s(i int) float64 {
	if i == 0 {
		return 0.0
	}
	r := rsamples[i]
	return -math.Sqrt(3) / 9.0 * (2.0*zEff4*r*r/243.0 - 2.0*zEff3*r/9.0 + 13.0*zEff2/9.0 - 2.0*zEff/r) * math.Exp(-zEff*r/3.0)
}

// Generates the 3x3 hamiltonian matrix in the {|1s〉|2s〉|3s〉} basis.
func hamiltonian(n, uEE func(int) float64) *mat.SymDense {
	calc := func(bra, ket, kinKet func(int) float64) float64 {
		return zEff3 * radialIntegrate(func(i int) float64 {
			braVal := bra(i)
			return (braVal*kinKet(i) + braVal*vKS(i, uEE, n)*ket(i))
		})
	}

	h11 := calc(r1s, r1s, kinR1s)
	h12 := calc(r1s, r2s, kinR2s)
	h13 := calc(r1s, r3s, kinR3s)

	h21 := calc(r2s, r1s, kinR1s)
	h22 := calc(r2s, r2s, kinR2s)
	h23 := calc(r2s, r3s, kinR3s)

	h31 := calc(r3s, r1s, kinR1s)
	h32 := calc(r3s, r2s, kinR2s)
	h33 := calc(r3s, r3s, kinR3s)

	return mat.NewSymDense(3, []float64{
		h11, h12, h13,
		h21, h22, h23,
		h31, h32, h33,
	})
}

func makeUEEFunc(uEE []float64) func(int) float64 {
	return func(i int) float64 { return uEE[i] }
}

func main() {
	// Initial density function guess: Sakurai eq (7.104).
	n := func(i int) float64 {
		r := rsamples[i]
		return zEff3 / (4.0 * math.Pi) * math.Exp(-zEff*r)
	}

	// Generates the .m file for plotting with octave.
	filename := filepath.Join(os.TempDir(), "sakurai_dft.m")
	fmt.Println(filename)
	f, err := os.Create(filename)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	type entry struct {
		uExt float64
		tKS  float64
		uEE  float64
		uXC  float64
		ev   [3]float64
	}

	const iterations = 4
	entries := make([]entry, iterations+1)
	for k := 0; k <= iterations; k++ {
		var nSlice []float64
		rlimit := 6.0
		for ii := 0; rsamples[ii] <= rlimit; ii++ {
			r := rsamples[ii]
			nSlice = append(nSlice, r*r*n(ii))
		}
		fmt.Fprintf(f, "r%v=%v;\n", k, rsamples[:len(nSlice)])
		fmt.Fprintf(f, "n%v=%v;\n", k, nSlice)
		fmt.Fprintf(f, "plot(r%v,n%v,';r^2n^{(%v)};');hold on;\n", k, k, k)

		uEE := makeUEEFunc(uEESlice(n))

		// Uext: see eq (7.100).
		entries[k].uExt = 4.0 * math.Pi * radialIntegrate(func(i int) float64 {
			return vExt(i) * n(i)
		})

		// Uee: see eq (7.94).
		entries[k].uEE = 2.0 * math.Pi * radialIntegrate(func(i int) float64 {
			return n(i) * uEE(i)
		})

		// Uxc: see eq (7.101a).
		entries[k].uXC = 4.0 * math.Pi * radialIntegrate(func(i int) float64 {
			nVal := n(i)
			return nVal * epsilonXC(nVal)
		})

		h := hamiltonian(n, uEE)
		var eig mat.EigenSym
		ok := eig.Factorize(h, true)
		if !ok {
			panic("eigen decomposition failed")
		}
		var ev mat.Dense
		eig.VectorsTo(&ev)

		if k != iterations {
			entries[k+1].ev[0] = ev.At(0, 0)
			entries[k+1].ev[1] = ev.At(1, 0)
			entries[k+1].ev[2] = ev.At(2, 0)

			phiBra := func(i int) float64 {
				return ev.At(0, 0)*r1s(i) + ev.At(1, 0)*r2s(i) + ev.At(2, 0)*r3s(i)
			}
			kinPhiKet := func(i int) float64 {
				return ev.At(0, 0)*kinR1s(i) + ev.At(1, 0)*kinR2s(i) + ev.At(2, 0)*kinR3s(i)
			}
			// Tks: see eq (7.92).
			entries[k+1].tKS = 2.0 * zEff3 * radialIntegrate(func(i int) float64 {
				return phiBra(i) * kinPhiKet(i)
			})
		}

		// New iteration of density function, see eq (7.109).
		n = func(i int) float64 {
			t := ev.At(0, 0)*r1s(i) + ev.At(1, 0)*r2s(i) + ev.At(2, 0)*r3s(i)
			return t * t * zEff3 / (2.0 * math.Pi)
		}
	}

	fmt.Println("Iter    Uext        Tks        Uee        Uxc         E            Φ")
	for k := 1; k <= iterations; k++ {
		energy := (entries[k].uExt + entries[k].tKS + entries[k].uEE + entries[k].uXC) * hatree
		fmt.Printf(
			"%v       %1.5f    %1.5f    %1.5f    %1.5f    %1.5f    [%1.5f %1.5f %1.5f]\n",
			k, entries[k].uExt, entries[k].tKS, entries[k].uEE, entries[k].uXC, energy,
			entries[k].ev[0],
			entries[k].ev[1],
			entries[k].ev[2],
		)
	}

	fmt.Fprintf(f, "xlabel('r'); ylabel('r^2n(r)'); set(gca,'fontsize',16);")
}
