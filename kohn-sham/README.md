# numerically calculates helium atom ground state based on density function theory and kohn-sham equation

This program carries out all the detailed calculation for helium atom based on DFT and kohn-sham equation following the recipe given in section 7.6.5. Table 7.1 and Figure 7.7 are reproduced.

## Example - 4 iterations
```
./kohn-sham (master) ▶ go run main.go
/var/folders/_0/2d8v_l8x5r947l5f35hdx0yw0000gq/T/sakurai_dft.m
Iter    Uext        Tks        Uee        Uxc         E            Φ
1       -6.90899    2.98841    2.17556    -1.13173    -78.28018    [-0.99851 -0.04839 -0.02523]
2       -6.48463    2.63623    1.99610    -1.04550    -78.85316    [-0.99712 0.07204 0.02352]
3       -6.56551    2.69842    2.03065    -1.06192    -78.86826    [-0.99852 0.05225 0.01540]
4       -6.54931    2.68583    2.02373    -1.05863    -78.86890    [-0.99827 0.05632 0.01701]
```

Then use `octave` to plot with the generated .m file
```
./kohn-sham (master) ▶ octave --persist /var/folders/_0/2d8v_l8x5r947l5f35hdx0yw0000gq/T/sakurai_dft.m
```
<img width="1535" alt="Screenshot 2023-03-22 at 09 49 50" src="https://user-images.githubusercontent.com/107862003/226781000-1472270b-0f32-4cfa-9aff-4d416c2f7924.png">
