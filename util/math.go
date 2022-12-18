package util

import (
	"fmt"
	"math/big"
)

// Creates a brand new blank rational.
func BlankRat() *big.Rat { return big.NewRat(0, 1) }

// Creates a brand new blank int.
func BlankInt() *big.Int { return big.NewInt(0) }

// Returns factorial of z.
func Factorial(z int) *big.Int {
	f := big.NewInt(1)
	for k := 2; k <= z; k++ {
		f.Mul(f, big.NewInt(int64(k)))
	}
	return f
}

// Takes the square root of products of all z's, returns (a, b) for the result a√b.
func SqrtProducts(vals []int) (*big.Int, *big.Int) {
	factors := make(map[int]int)
	for _, z := range vals {
		if z < 0 {
			panic(fmt.Sprintf("SqrtProducts called with negative value %v", z))
		}
		terms := Factorize(z)
		for _, t := range terms {
			factors[t.Prime] += t.Power
		}
	}
	in, out := big.NewInt(1), big.NewInt(1)
	for prime, power := range factors {
		exp := BlankInt().Exp(big.NewInt(int64(prime)), big.NewInt(int64(power/2)), nil)
		out.Mul(out, exp)
		if power%2 == 1 {
			in.Mul(in, big.NewInt(int64(prime)))
		}
	}
	return in, out
}

// Term represents a factorization term
type Term struct {
	Prime int
	Power int
}

// Factorize factorizes z into prime factors. Just brute-force calculation, no optimizations whatsoever.
func Factorize(z int) []Term {
	var terms []Term
	var ct *Term
	q := (z + 1) / 2
	p := 2
	for p <= q {
		if z%p == 0 {
			if ct == nil {
				terms = append(terms, Term{Prime: p})
				ct = &terms[len(terms)-1]
			}
			ct.Power++
			z /= p
		} else {
			ct = nil
			p++
		}
	}
	if len(terms) == 0 {
		terms = []Term{{Prime: z, Power: 1}}
	}
	return terms
}
