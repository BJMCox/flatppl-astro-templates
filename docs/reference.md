# Function reference

All arrays are real; matching vectors use the same positional grid. Physical
densities have declared units per coordinate or solid angle. Positive quadrature
weights carry the coordinate measure; already integrated bins must not receive
them twice. Domains below are caller obligations; these arithmetic functions
do not validate or clip inputs.

Choose numerical units that keep both values and fitted-parameter derivatives
within the executing precision. Raw SI luminosity-distance derivatives can be
subnormal in float32 even when the flux is representable. For
$L_\lambda=10^{30}a$ W/m and $D_L=10^{20}d$ m, call the wavelength helper
with numerical values `a, z, d` in those units, then multiply its result by
`1e-10` to obtain W/m²/m. This preserves gradients in the dimensionless fit
parameters, including zero strength; see
[the numerical-units query](../tests/runtime-distance-scaled.flatppl).
The physical distance law remains the same. Raw SI values are covered by
[the distance-range query](../tests/runtime-distance.flatppl); use float64
when SI adjoints fall outside the float32 range.

| Function | Result |
|---|---|
| [`redshift_wavelength`](#redshift_wavelength) | Observed to rest wavelength |
| [`redshift_flux_lambda`](#redshift_flux_lambda) | Rest luminosity per wavelength to observed flux per wavelength |
| [`redshift_frequency`](#redshift_frequency) | Observed to rest frequency |
| [`redshift_flux_nu`](#redshift_flux_nu) | Rest luminosity per frequency to observed flux per frequency |
| [`dust_transmission`](#dust_transmission) | Magnitude extinction to dimensionless transmission |
| [`photon_band_counts`](#photon_band_counts) | Expected photon counts in one passband |
| [`energy_band_flux`](#energy_band_flux) | Integrated transmitted energy flux |
| [`normalize_sky`](#normalize_sky) | Sky density per steradian |

<a id="redshift_wavelength"></a>
## `redshift_wavelength(wavelength, redshift)`

`wavelength` is a positive observed-frame scalar or vector; `redshift` is a
scalar greater than -1. Returns rest-frame coordinates in the same units and shape:

```math
\lambda_{\rm rest}=\frac{\lambda_{\rm obs}}{1+z}.
```

Evaluate the rest-frame template at these coordinates before converting its density.

<a id="redshift_flux_lambda"></a>
## `redshift_flux_lambda(rest_luminosity_density, redshift, luminosity_distance)`

`rest_luminosity_density` is nonnegative scalar or vector $L_\lambda$, already
evaluated at $\lambda_{\rm obs}/(1+z)$. `redshift` is scalar $z>-1$;
`luminosity_distance` is scalar $D_L>0$. Returns the same shape:

```math
F_\lambda(\lambda_{\rm obs})=\frac{L_\lambda(\lambda_{\rm obs}/(1+z))}{4\pi D_L^2(1+z)}.
```

For $L_\lambda$ in W/m and $D_L$ in m, the result is W/m²/m. Distance and
redshift are explicit; their physical relation is the caller's cosmology.
Since $d\lambda_{\rm obs}=(1+z)d\lambda_{\rm rest}$, the corresponding-domain
integral is $F_{\rm bol}=L_{\rm bol}/(4\pi D_L^2)$. Do not shift coordinates twice.

<a id="redshift_frequency"></a>
## `redshift_frequency(frequency, redshift)`

`frequency` is a positive observed-frame scalar or vector; `redshift` is scalar
$z>-1$. Returns rest-frame coordinates in the same units and shape:

```math
\nu_{\rm rest}=(1+z)\nu_{\rm obs}.
```

<a id="redshift_flux_nu"></a>
## `redshift_flux_nu(rest_luminosity_density, redshift, luminosity_distance)`

`rest_luminosity_density` is nonnegative scalar or vector $L_\nu$, already
evaluated at $(1+z)\nu_{\rm obs}$. `redshift` is scalar $z>-1$;
`luminosity_distance` is scalar $D_L>0$. Returns the same shape:

```math
F_\nu(\nu_{\rm obs})=\frac{(1+z)L_\nu((1+z)\nu_{\rm obs})}{4\pi D_L^2}.
```

For W/Hz and m inputs the result is W/m²/Hz. Here
$d\nu_{\rm obs}=d\nu_{\rm rest}/(1+z)$, so its integral gives the same
bolometric distance law. Wavelength and frequency densities require different
Jacobians; when converting them use $F_\nu=F_\lambda\lambda^2/c$ with compatible units.

<a id="dust_transmission"></a>
## `dust_transmission(extinction_magnitudes)`

`extinction_magnitudes` is nonnegative scalar or array $A$ in magnitudes.
Returns dimensionless transmission of the same shape, in $(0,1]$ mathematically:

```math
T=10^{-0.4A}.
```

The caller evaluates its extinction curve in the declared rest or observed
frame and multiplies the relevant density by $T$. No extinction curve or dust
emission is inferred. Very large extinction can underflow numerically to zero.

<a id="photon_band_counts"></a>
## `photon_band_counts(wavelengths, flux_lambda, effective_area, quadrature_weights, exposure, planck_constant, speed_of_light)`

`wavelengths`, `flux_lambda`, `effective_area`, and `quadrature_weights` are
aligned vectors `[N]`: positive observed wavelengths, nonnegative energy flux
per wavelength after attenuation, nonnegative effective collecting areas
including throughput, and positive wavelength quadrature widths. `exposure`
is a nonnegative scalar duration; `planck_constant` and `speed_of_light` are
positive scalar $h,c$ in compatible units. Returns a scalar expected count:

```math
C=t\sum_i A_{{\rm eff},i}F_{\lambda,i}q_i\frac{\lambda_i}{hc}.
```

For wavelengths and widths in m, W/m²/m flux, m² area, seconds, J s, and m/s,
the result is dimensionless photon counts. Throughput, area, exposure, and bin
widths each enter once. A passband with zero throughput yields zero counts.
This is a mean; the caller specifies a sampling law.

<a id="energy_band_flux"></a>
## `energy_band_flux(flux_lambda, transmission, quadrature_weights)`

Aligned vectors `[N]` contain nonnegative energy flux per wavelength,
dimensionless transmission in `[0,1]`, and positive wavelength quadrature widths.
Returns a scalar integrated transmitted energy flux:

```math
F_{\rm band}=\sum_i F_{\lambda,i}T_iq_i.
```

W/m²/m and m yield W/m². This integral has no photon weighting or normalization
to a band average. A previously attenuated flux must not receive the same dust
transmission again. Calibrated averages and photometric zero points need their
own explicitly stated convention.

<a id="normalize_sky"></a>
## `normalize_sky(intensity, solid_angles)`

Aligned vectors `[N]` contain nonnegative intensity $I_i$ per steradian and
positive pixel solid angles $\Omega_i$ in steradians. The weighted sum must be
positive. Returns a vector density per steradian on the supplied sky domain:

```math
S_i=\frac{I_i}{\sum_j I_j\Omega_j},\qquad \sum_i S_i\Omega_i=1.
```

Multiplying by $\Omega_i$ gives pixel fractions. Unequal-area pixels require
their individual areas. The caller states frame, projection, and mask; a sky
subset is normalized on that subset. Zero-intensity pixels stay zero. This is
a shape normalization and removes the original integrated strength.

## Composition and mathematical boundaries

Use [flatppl-stats-templates](https://github.com/BJMCox/flatppl-stats-templates)
for template combination, weighted projection, interpolation, and observation
models. Templates have bins as rows and components as columns `[N,K]`;
responses are `[M,N]`. Physical strengths need not sum to one. Fixed responses
can be folded once; fitted attenuation or calibration requires recomputation.
Losses outside a selected response domain stay explicit.

Velocity broadening is not implemented in 0.1. A flux-preserving convolution
on $u=\log\lambda$ would act on $G(u)=\lambda F_\lambda(\lambda)$, with
$G'(u)=\int k(u-v)G(v)dv$ and $\int k(u)du=1$, then return $F'_\lambda=G'/\lambda$.
Convolving $F_\lambda$ directly on this grid does not in general preserve
$\int F_\lambda d\lambda$. A discrete operator must state quadrature weights,
kernel boundaries, and any losses; renormalizing away those losses changes
the model. No finite-grid broadening convention is silently selected here.

## Synthetic example ABI

[spectrum-photometry.flatppl](../examples/spectrum-photometry.flatppl) accepts
one nonnegative strength vector `[2]`. It returns attenuated flux `[3]` in
W/m²/m, source counts `[2]`, total counts `[2]`, and the scalar Poisson log
likelihood for the inline synthetic observations. Constants state the observed
wavelength grid, widths, dust magnitudes, effective areas, seconds, $h$, and $c$.
The data constrain combinations of physical strengths; no prior or population
synthesis model is supplied.
For strengths `[1.2, 0.7]`, total means are approximately
`[71.4190276832, 89.9814951027]`, and the log likelihood for `[70, 90]` is
`-6.22833948179`, including Poisson factorial constants.
