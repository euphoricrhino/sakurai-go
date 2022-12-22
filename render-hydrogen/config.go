package main

import (
	"encoding/json"
	"fmt"
	"image/color"
	"image/png"
	"math"
	"os"
)

// Configures the hydrogen renderer.
type config struct {
	// Theta angle of the camera, in degrees.
	CameraTheta float64 `json:"cameraTheta"`
	// Phi angle of the camera, in degrees,
	CameraPhi float64 `json:"cameraPhi"`
	// Image size in pixels.
	ImageSize int `json:"imageSize"`
	// Field-of-view size in units of bohr radius.
	FOVSize float64 `json:"fovSize"`
	// Distance between sampling layers.
	LayerDist float64 `json:"layerDist"`
	// Sampling layers to compute, must be positive odd number, the middle of which represents the perpendicular plane of sight containing origin.
	Layers int `json:"layers"`

	Concurrency int `json:"concurrency"`
	// Precision in binary digits for the float.
	FloatPrec uint `json:"floatPrec"`

	// .PNG file for the heatmap.
	HeatmapFile string  `json:"heatmapFile"`
	OutputFile  string  `json:"outputFile"`
	Exposure    float32 `json:"exposure"`
	// Quantum numbers.
	N int `json:"n"`
	L int `json:"l"`
	M int `json:"m"`

	heatmap []color.Color
}

const degToRad = math.Pi / 180.0

func parseConfigOrDie(filename string) *config {
	cfg := &config{}
	data, err := os.ReadFile(filename)
	if err != nil {
		panic(err)
	}
	if err := json.Unmarshal(data, cfg); err != nil {
		panic(fmt.Sprintf("failed to unmarshal config: %v", err))
	}
	// Convert angles into radians.
	cfg.CameraTheta *= degToRad
	cfg.CameraPhi *= degToRad

	if cfg.ImageSize <= 1 {
		panic(fmt.Sprintf("invalid imageSize: %v", cfg.ImageSize))
	}
	if cfg.FOVSize <= 0 {
		panic(fmt.Sprintf("invalid fovSize: %v", cfg.FOVSize))
	}
	if cfg.LayerDist <= 0 {
		panic(fmt.Sprintf("invalid layerDist: %v", cfg.LayerDist))
	}
	if cfg.Layers <= 0 || cfg.Layers%2 == 0 {
		panic(fmt.Sprintf("invalid layers: %v", cfg.Layers))
	}
	if cfg.Concurrency <= 0 {
		panic(fmt.Sprintf("invalid concurrency: %v", cfg.Concurrency))
	}
	if cfg.Exposure <= 0 {
		panic(fmt.Sprintf("invalid exposure: %v", cfg.Exposure))
	}

	// Validate quantum numbers.
	if cfg.N <= 0 {
		panic(fmt.Sprintf("invalid n: %v", cfg.N))
	}
	if cfg.L < 0 || cfg.L >= cfg.N {
		panic(fmt.Sprintf("invalid l: %v", cfg.L))
	}
	if cfg.M < 0 || cfg.M > cfg.L {
		panic(fmt.Sprintf("invalid m: %v", cfg.M))
	}

	// Load heatmap file.
	f, err := os.Open(cfg.HeatmapFile)
	if err != nil {
		panic(fmt.Sprintf("failed to open heatmap file: %v", err))
	}
	defer f.Close()
	hm, err := png.Decode(f)
	if err != nil {
		panic(fmt.Sprintf("failed to decode heatmap as PNG: %v", err))
	}
	rect := hm.Bounds()
	width := rect.Max.X - rect.Min.X
	cfg.heatmap = make([]color.Color, width)
	for i := 0; i < width; i++ {
		cfg.heatmap[i] = hm.At(i+rect.Min.X, rect.Min.Y)
	}
	return cfg
}
