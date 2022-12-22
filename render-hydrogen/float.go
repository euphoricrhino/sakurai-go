package main

import (
	"math/big"
	"sync"
)

var (
	// Precision of all big.Float computation.
	floatPrec uint = 2000
	once      sync.Once
)

// Will be set only once after config is loaded.
func setPrecOnce(prec uint) {
	once.Do(func() {
		floatPrec = prec
	})
}

func blankFloat() *big.Float { return big.NewFloat(0).SetPrec(floatPrec) }

func newFromFloat64(val float64) *big.Float {
	return big.NewFloat(val).SetPrec(floatPrec)
}

func newFromInt(val int) *big.Float {
	return blankFloat().SetInt64(int64(val)).SetPrec(floatPrec)
}

func newFromRat(n, d int) *big.Float {
	return blankFloat().SetRat(big.NewRat(int64(n), int64(d))).SetPrec(floatPrec)
}
