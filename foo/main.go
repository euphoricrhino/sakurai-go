package main

import (
	"flag"
	"fmt"
	"math/big"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"text/template"
)

func main() {
	var bstr string

	flag.StringVar(&bstr, "b", "", "b")

	flag.Parse()

	strs := strings.Split(bstr, ",")
	var bs []int
	for _, p := range strs {
		v, err := strconv.Atoi(p)
		if err != nil {
			panic(err)
		}
		bs = append(bs, v)
	}

	str := "\\begin{aligned} "
	for _, b := range bs {
		str += formula(b)
		str += "\\\\"
	}
	str += "\\end{aligned}"
	renderMath(str, "2023.html")
}

const target = 2023

var (
	one = big.NewInt(1)
)

func formula(b int) string {
	if b < (target*target)/4 {
		panic(fmt.Sprintf("too small: %v", b))
	}
	p := b - target*target
	q := -target * b

	s := big.NewRat(int64(-q), 2)
	qq := big.NewInt(int64(q))
	qq.Mul(qq, qq)
	t := big.NewRat(1, 1).SetFrac(qq, big.NewInt(4))
	ppp := big.NewInt(int64(p))
	ppp.Mul(ppp, big.NewInt(int64(p)))
	ppp.Mul(ppp, big.NewInt(int64(p)))
	t.Add(t, big.NewRat(1, 1).SetFrac(ppp, big.NewInt(27)))

	sstr := ""
	if s.Denom().Cmp(one) == 0 {
		sstr = s.Num().String()
	} else {
		sstr = fmt.Sprintf("\\frac{%v}{%v}", s.Num(), s.Denom())
	}

	tstr := ""
	if t.Denom().Cmp(one) == 0 {
		tstr = t.Num().String()
	} else {
		tstr = fmt.Sprintf("\\frac{%v}{%v}", t.Num(), t.Denom())
	}
	return fmt.Sprintf("%v&=\\sqrt[3]{%v+\\sqrt{%v}}+\\sqrt[3]{%v-\\sqrt{%v}}", target, sstr, tstr, sstr, tstr)
}

const tmplStr = `<!DOCTYPE html>
<html>
<head>
<script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
<script type="text/javascript" id="MathJax-script" async
  src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js">
</script>
</head>
<body>
$${{ . }}$$
</body>
</html>`

var tmpl *template.Template

func init() {
	tmpl = template.Must(template.New("root").Parse(tmplStr))
}

func renderMath(latexStr string, toFile string) {
	filename := filepath.Join(os.TempDir(), toFile)
	f, err := os.Create(filename)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	if err := tmpl.Execute(f, latexStr); err != nil {
		panic(err)
	}

	fmt.Println(filename)
}
