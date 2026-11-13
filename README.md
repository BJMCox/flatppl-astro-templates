# flatppl-astro-templates

[Function reference](docs/reference.md)

FlatPPL 0.1 operators for redshifted spectral densities, dust attenuation,
photon and energy passbands, and solid-angle sky normalization.

Run the [posterior example](https://live.flatppl.org/#model=https%3A%2F%2Fraw.githubusercontent.com%2Fflatppl-community%2Fflatppl-astro-templates%2Fmain%2Fexamples%2Fplayground.flatppl&target=Pi_post) in the live playground. It loads its library modules from GitHub.

```flatppl
astro = load_module("main.flatppl")
rest_wavelength = astro.redshift_wavelength([4e-7, 5e-7], 0.2)
transmission = astro.dust_transmission([0.2, 0.1])
```

The [function reference](docs/reference.md) defines units, axes, Jacobians,
and domains. The [synthetic spectrum/photometry example](examples/spectrum-photometry.flatppl)
combines supplied physical spectra and projects them into photon counts using
[flatppl-stats-templates](https://github.com/flatppl-community/flatppl-stats-templates).
Its shared dependency is pinned in the example's `load_module` URL.

Every physical input is explicit, including luminosity distance, Planck's
constant, light speed, quadrature widths, and effective area. Use compatible
units. The library supplies no cosmology, stellar-population generator, dust
emission, implicit prior, or density floor. Response losses remain losses.

For local validation with an existing Node and NumPy/SciPy environment:

```sh
FLATPPL_JS_DIR=/path/to/flatppl-js python tests/check.py /path/to/flatppl-stats-templates
```

The tests use the JS module loader and independent analytic/NumPy/SciPy
oracles. [tests/runtime.flatppl](tests/runtime.flatppl) and the example have
explicit `inputs` and `outputs` for StableHLO compilation.
