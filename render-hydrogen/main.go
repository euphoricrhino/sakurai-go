package main

import (
	"flag"
	"fmt"
	"image"
	"image/color"
	"image/png"
	"math"
	"math/big"
	"os"
	"sync"
)

func main() {
	var configFile string
	flag.StringVar(&configFile, "config", "", "config file")
	flag.Parse()

	cfg := parseConfigOrDie(configFile)

	setPrecOnce(cfg.FloatPrec)

	scr := newScreen(cfg)
	eval := newEvaluator(cfg)

	// Data with all layers added together.
	data := make([]*big.Float, cfg.ImageSize*cfg.ImageSize)
	var wg sync.WaitGroup
	wg.Add(cfg.Concurrency + 1)

	ch := make(chan struct{}, cfg.Concurrency)

	// Counting goroutine, will display percentage progress.
	go func() {
		defer wg.Done()
		mark := 0
		cnt := 0
		total := len(data) * cfg.Layers
		for cnt < total {
			<-ch
			cnt++
			progress := int(1000.0 * float64(cnt) / float64(total))
			if progress > mark {
				mark = progress
				fmt.Printf("\r       \r%03.1f%%", float64(progress)/10.0)
			}
		}
		fmt.Printf("\r       \rdone\n")
	}()

	// One goroutine per worker, with all workers partitioning the i-index space.
	for w := 0; w < cfg.Concurrency; w++ {
		go func(shard int) {
			defer wg.Done()
			for i := 0; i < cfg.ImageSize; i++ {
				if i%cfg.Concurrency != shard {
					continue
				}
				for j := 0; j < cfg.ImageSize; j++ {
					// Recall cfg.Layers must be odd number
					halfLayers := (cfg.Layers - 1) / 2
					pixel := blankFloat()
					for k := -halfLayers; k <= halfLayers; k++ {
						r := scr.gridToWorld(i, j, k)
						p := eval.probDensity(r[0], r[1], r[2])
						pixel.Add(pixel, p)
						ch <- struct{}{}
					}
					// Safe to write since no other worker goroutine will touch the same (i,j).
					data[j*cfg.ImageSize+i] = pixel
				}
			}
		}(w)
	}

	wg.Wait()

	// For normalization calculation.
	max := blankFloat()
	for i := 0; i < len(data); i++ {
		if max.Cmp(data[i]) < 0 {
			max.Set(data[i])
		}
	}

	img := image.NewRGBA(image.Rect(0, 0, cfg.ImageSize, cfg.ImageSize))
	for i := 0; i < cfg.ImageSize; i++ {
		for j := 0; j < cfg.ImageSize; j++ {
			pixel := data[j*cfg.ImageSize+i]
			normalized, _ := blankFloat().Quo(pixel, max).Float64()
			// Adjust for exposure for best visual contrast.
			val := math.Pow(normalized, 1.0/float64(cfg.Exposure))
			heatmapPos := int(val * float64(len(cfg.heatmap)-1))
			r, g, b, a := cfg.heatmap[heatmapPos].RGBA()
			img.SetRGBA64(i, j, color.RGBA64{R: uint16(r), G: uint16(g), B: uint16(b), A: uint16(a)})
		}
	}
	out, err := os.Create(cfg.OutputFile)
	if err != nil {
		panic(err)
	}
	defer out.Close()
	if err := png.Encode(out, img); err != nil {
		panic(err)
	}
}
