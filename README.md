# flatppl-astro-templates

[Function reference](docs/reference.md) · [Live example](https://live.flatppl.org/#model=https%3A%2F%2Fraw.githubusercontent.com%2Fflatppl-community%2Fflatppl-astro-templates%2Fmain%2Fexamples%2Fplayground.flatppl&target=Pi_post)

Redshifted spectral densities, dust attenuation, photon and energy passbands, and solid-angle sky normalization for FlatPPL. Library version **0.1**.

```flatppl
astro = load_module("main.flatppl")
rest_wavelength = astro.redshift_wavelength([4e-7, 5e-7], 0.2)
transmission = astro.dust_transmission([0.2, 0.1])
```

Physical constants, luminosity distance, quadrature widths, and effective area are explicit inputs. Units must agree. Response losses remain losses.

The [spectrum and photometry example](examples/spectrum-photometry.flatppl) projects supplied spectra into photon counts.

See the reference for definitions, domains, numerical limits, and validation commands.
