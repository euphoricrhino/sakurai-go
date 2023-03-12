# plots resonance scattering of partial wave for constant spherical well potential
The go program will generate a .m file which can be passed to octave for plotting.

## Example - resonance scattering
```
./resonance-scattering (master) ▶ octave --persist `go run *.go --l=3 --x-start=1.3 --x-step=.001 -x-limit=1.5 --gamma=5.5`
```
will reproduce Figure 6.15
<img width="1313" alt="Screenshot 2023-03-12 at 21 20 04" src="https://user-images.githubusercontent.com/107862003/224547298-69d330f8-34bd-40de-adbe-d9993d872ad6.png">
Some other interesting combinations:
```
./resonance-scattering (master) ▶ octave --persist `go run *.go --l=23 --x-start=18 --x-step=.001 -x-limit=19 --gamma=20.5 --prec=2000`
```
<img width="1011" alt="Screenshot 2023-03-12 at 21 20 52" src="https://user-images.githubusercontent.com/107862003/224547342-f9705703-2f95-4de7-bc8b-dbe70d2e3a74.png">

```
./resonance-scattering (master) ▶ octave --persist `go run *.go --l=35 --x-start=10000 --x-step=1 -x-limit=100000 --gamma=909.5 --prec=2000`
```

<img width="1010" alt="Screenshot 2023-03-12 at 21 22 30" src="https://user-images.githubusercontent.com/107862003/224547459-8fc8fe52-0314-4bd9-9f99-0fd52864a0ad.png">
