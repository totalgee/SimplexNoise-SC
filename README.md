# SimplexNoise Quark

This SuperCollider [Quark](http://doc.sccode.org/Guides/UsingQuarks.html)
provides a set of helpers to generate
[Simplex Noise](https://en.wikipedia.org/wiki/Simplex_noise)
in 2D and 3D space.

The `Simplex.noise2` and `Simplex.noise3` functions are based on
the example Java code (in the public domain) by Stefan Gustavson,
with optimizations by Peter Eastman. You can read
[the paper](http://staffwww.itn.liu.se/~stegu/simplexnoise/simplexnoise.pdf)
or see the [original code](http://weber.itn.liu.se/~stegu/simplexnoise/SimplexNoise.java).

Besides the underlying Simplex implementation, there are extra
helpers (`Simplex.fBm2` and `Simplex.fBm3`) for
[fractal/fractional Brownian motion](https://en.wikipedia.org/wiki/Fractional_Brownian_motion)
(i.e. summed octaves of noise) in 2D and 3D, as well as
1D periodic noise (`Simplex.periodic`).

## Getting started
To install SimplexNoise as a Quark in SuperCollider, run the
following line:

```supercollider
Quarks.install("https://github.com/totalgee/SimplexNoise-SC");
```
