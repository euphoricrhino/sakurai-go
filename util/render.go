package util

import (
	"fmt"
	"html/template"
	"os"
	"path/filepath"
)

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

func RenderMath(latexStr string, toFile string) {
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
