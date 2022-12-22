package main

import (
	"math"
	"math/big"
)

// Calculates screen to world geometry.
type screen struct {
	// Normalized "in" vector from center-of-screen to the origin, in world coordinates.
	in [3]*big.Float
	// Normalized "up" vector of screen, in world coordinates.
	up [3]*big.Float
	// Normalized "right" vector of screen, in world coordinates.
	right [3]*big.Float
	// Middle layer's northwest corner pixel, in world coordinate.
	nw [3]*big.Float
	// Steps in world measurement to the right/up/in direction respectively.
	step [3]*big.Float
}

func newScreen(cfg *config) *screen {
	ct, st := math.Cos(cfg.CameraTheta), math.Sin(cfg.CameraTheta)
	cp, sp := math.Cos(cfg.CameraPhi), math.Sin(cfg.CameraPhi)

	pixelStep := newFromFloat64(cfg.FOVSize / float64(cfg.ImageSize-1))
	scr := &screen{
		in: [3]*big.Float{
			newFromFloat64(-st * cp),
			newFromFloat64(-st * sp),
			newFromFloat64(-ct),
		},
		up: [3]*big.Float{
			newFromFloat64(-ct * cp),
			newFromFloat64(-ct * sp),
			newFromFloat64(st),
		},
		right: [3]*big.Float{
			newFromFloat64(-sp),
			newFromFloat64(cp),
			blankFloat(),
		},
		step: [3]*big.Float{
			pixelStep,
			pixelStep,
			newFromFloat64(cfg.LayerDist),
		},
	}
	hs := newFromFloat64(cfg.FOVSize * 0.5)
	scr.nw[0] = blankFloat().Mul(hs, blankFloat().Sub(scr.up[0], scr.right[0]))
	scr.nw[1] = blankFloat().Mul(hs, blankFloat().Sub(scr.up[1], scr.right[1]))
	scr.nw[2] = blankFloat().Mul(hs, blankFloat().Sub(scr.up[2], scr.right[2]))

	return scr
}

// Computes world coordinates given grid coordinate (i,j,k), where (i,j) is the x-y pixel coordinate
// and k is the layer index (middle layer contains origin).
func (scr *screen) gridToWorld(i, j, k int) [3]*big.Float {
	ret := [3]*big.Float{}
	for n := 0; n < 3; n++ {
		ret[n] = blankFloat().Set(scr.nw[n])

		right := newFromInt(i)
		right.Mul(right, scr.step[0])
		right.Mul(right, scr.right[n])
		ret[n].Add(ret[n], right)

		up := newFromInt(j)
		up.Mul(up, scr.step[1])
		up.Mul(up, scr.up[n])
		ret[n].Sub(ret[n], up)

		in := newFromInt(k)
		in.Mul(in, scr.step[2])
		in.Mul(in, scr.in[n])
		ret[n].Add(ret[n], in)
	}
	return ret
}
