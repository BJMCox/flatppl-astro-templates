"""Independent physical invariants through the JS module loader; no downloads."""

import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
from numpy.testing import assert_allclose
from scipy.integrate import quad
from scipy.stats import poisson

ROOT = Path(__file__).resolve().parents[1]
if len(sys.argv) > 1:
    os.environ["FLATPPL_STATS_TEMPLATES_DIR"] = str(Path(sys.argv[1]).resolve())


def evaluate(request):
    result = subprocess.run(
        ["node", str(ROOT / "tests/evaluate.cjs")],
        input=json.dumps(request), text=True, capture_output=True, check=False,
    )
    if result.returncode:
        raise RuntimeError(result.stderr)
    return json.loads(result.stdout)


def literal(values):
    return json.dumps(np.asarray(values).tolist())


def runtime_oracle(x):
    """Independent NumPy expression for tests/runtime.flatppl, x shape [4]."""
    z, distance, extinction, amplitude = x
    wavelength = np.array([4e-7, 5e-7, 6e-7])
    luminosity = 1e-9 * amplitude * np.array([2.0, 3.0, 1.0])
    flux = luminosity / (4 * np.pi * distance**2 * (1 + z))
    flux_nu = luminosity * (1 + z) / (4 * np.pi * distance**2)
    transmission = np.exp(-extinction * np.log(10) / 2.5)
    photon_energy = 6.62607015e-34 * 299792458 / wavelength
    photon_count = 100 * np.dot(flux * transmission / photon_energy, np.array([0.8, 0.5, 0.7]) * 1e-7)
    energy_flux = np.dot(flux, np.array([0.8, 0.5, 0.7]) * 1e-7)
    intensity = np.array([amplitude, 2.0, 0.0])
    sky = intensity / np.dot(intensity, [0.1, 0.3, 0.6])
    return np.concatenate((wavelength / (1 + z), np.array([3e14, 4e14, 5e14]) * (1 + z), flux, flux_nu, [transmission, photon_count, energy_flux], sky))


def distance_oracle(x):
    luminosity_lambda, luminosity_nu, distance, z = np.asarray(x, dtype=np.float64)
    area = 4 * np.pi * distance**2
    return np.array([luminosity_lambda / area / (1 + z), luminosity_nu * (1 + z) / area])


def runtime_jacobian(x):
    z, distance, _, amplitude = x
    value = runtime_oracle(x)
    jacobian = np.zeros((18, 4))
    jacobian[:3, 0] = -value[:3] / (1 + z)
    jacobian[3:6, 0] = value[3:6] / (1 + z)
    jacobian[6:9, 0] = -value[6:9] / (1 + z)
    jacobian[9:12, 0] = value[9:12] / (1 + z)
    jacobian[6:12, 1] = -2 * value[6:12] / distance
    jacobian[6:12, 3] = runtime_oracle([z, distance, x[2], 1.0])[6:12]
    jacobian[12, 2] = -np.log(10) / 2.5 * value[12]
    for index in (13, 14):
        jacobian[index, 0] = -value[index] / (1 + z)
        jacobian[index, 1] = -2 * value[index] / distance
        jacobian[index, 3] = runtime_oracle([z, distance, x[2], 1.0])[index]
    jacobian[13, 2] = -np.log(10) / 2.5 * value[13]
    normalization = 0.1 * amplitude + 0.6
    jacobian[15:18, 3] = np.array([0.6, -0.2, 0]) / normalization**2
    return jacobian


def distance_jacobian(x):
    _, _, distance, z = x
    value = distance_oracle(x)
    jacobian = np.zeros((2, 4))
    jacobian[0, 0] = 1 / (4 * np.pi * distance**2 * (1 + z))
    jacobian[1, 1] = (1 + z) / (4 * np.pi * distance**2)
    jacobian[:, 2] = -2 * value / distance
    jacobian[:, 3] = [-value[0] / (1 + z), value[1] / (1 + z)]
    return jacobian


def scaled_distance_oracle(x):
    strength, distance, z = x
    coefficients = np.array([1e-10 / (1 + z), 1e-25 * (1 + z)]) / (4 * np.pi * distance**2)
    value = strength * coefficients
    jacobian = np.column_stack((coefficients, -2 * value / distance, [-value[0] / (1 + z), value[1] / (1 + z)]))
    return value, jacobian


def check():
    # Coordinate Jacobians preserve the energy integral on transformed domains.
    z, distance = 0.7, 3.0
    rest = np.linspace(1.0, 4.0, 31)
    luminosity = 2 + 3 * rest
    observed_lambda = rest * (1 + z)
    observed_nu = rest / (1 + z)
    flux_lambda, flux_nu, coord_lambda, coord_nu = evaluate([
        f"a.redshift_flux_lambda({literal(luminosity)}, {z}, {distance})",
        f"a.redshift_flux_nu({literal(luminosity)}, {z}, {distance})",
        f"a.redshift_wavelength({literal(observed_lambda)}, {z})",
        f"a.redshift_frequency({literal(observed_nu)}, {z})",
    ])
    bolometric = quad(lambda v: 2 + 3 * v, 1, 4)[0] / (4 * np.pi * distance**2)
    assert_allclose(np.trapezoid(flux_lambda, observed_lambda), bolometric, rtol=2e-15, atol=0)
    assert_allclose(np.trapezoid(flux_nu, observed_nu), bolometric, rtol=2e-15, atol=0)
    assert_allclose([coord_lambda, coord_nu], [rest, rest], rtol=2e-15, atol=0)

    # The two densities describe the same physical spectrum under nu=c/lambda.
    lam, light = np.array([4e-7, 7e-7]), 299792458.0
    rest_lam = lam / (1 + z)
    l_lambda = np.array([2.0, 3.0])
    l_nu = l_lambda * rest_lam**2 / light
    f_lambda, f_nu = evaluate([
        f"a.redshift_flux_lambda({literal(l_lambda)}, {z}, {distance})",
        f"a.redshift_flux_nu({literal(l_nu)}, {z}, {distance})",
    ])
    assert_allclose(f_nu, np.array(f_lambda) * lam**2 / light, rtol=2e-15, atol=0)
    distance_query = (ROOT / "tests/runtime-distance.flatppl").read_text()
    distance_point = [1e30, 1e15, 1e20, 0.0]
    source = distance_query.replace("x = elementof(cartpow(reals, 4))", f"x = {literal(distance_point)}")
    densities = evaluate({"source": source, "bindings": ["result"]})[0]
    assert_allclose(densities, distance_oracle(distance_point), rtol=1e-14, atol=0)
    scaled_query = (ROOT / "tests/runtime-distance-scaled.flatppl").read_text()
    for point in ([1.0, 1.0, 0.0], [0.0, 1.0, 0.0]):
        source = scaled_query.replace("x = elementof(cartpow(reals, 3))", f"x = {literal(point)}")
        value = evaluate({"source": source, "bindings": ["result"]})[0]
        assert_allclose(value, scaled_distance_oracle(point)[0], rtol=2e-15, atol=0)

    # Photon and energy passbands differ by each photon's energy, not a band average.
    wavelengths = np.array([4e-7, 7e-7])
    flux, area, widths = np.array([2e-9, 3e-9]), np.array([0.3, 0.8]), np.array([1e-7, 2e-7])
    photon, energy, dust, zero_band = evaluate([
        f"a.photon_band_counts({literal(wavelengths)}, {literal(flux)}, {literal(area)}, {literal(widths)}, 120.0, 6.62607015e-34, {light})",
        f"a.energy_band_flux({literal(flux)}, {literal(area)}, {literal(widths)})",
        "a.dust_transmission([0.0, 1.0, 2.5])",
        "a.photon_band_counts([4e-7], [2e-9], [0.0], [1e-7], 120.0, 6.62607015e-34, 299792458.0)",
    ])
    photon_energy = 6.62607015e-34 * light / wavelengths
    assert_allclose(photon, 120 * np.dot(flux * widths, area / photon_energy), rtol=2e-15, atol=0)
    assert_allclose(energy, np.dot(flux * widths, area), rtol=2e-15, atol=0)
    assert_allclose(dust, np.exp(-np.array([0.0, 1.0, 2.5]) * np.log(10) / 2.5), rtol=2e-15, atol=0)
    assert zero_band == 0

    # Unequal-area sky pixels normalize physical area, preserve zeros and scale invariance.
    sky, scaled = evaluate([
        "a.normalize_sky([1.0, 2.0, 0.0], [0.1, 0.3, 0.6])",
        "a.normalize_sky([9.0, 18.0, 0.0], [0.1, 0.3, 0.6])",
    ])
    assert_allclose(sky, np.array([1, 2, 0]) / 0.7, rtol=2e-15, atol=0)
    assert_allclose(np.dot(sky, [0.1, 0.3, 0.6]), 1, rtol=2e-15, atol=0)
    assert_allclose(scaled, sky, rtol=2e-15, atol=0)
    assert sky[2] == 0

    # Runtime query and physical example are tested at concrete, independent points.
    runtime = (ROOT / "tests/runtime.flatppl").read_text()
    for point in ([0.7, 3.0, 0.2, 1.5], [0.0, 2.0, 0.0, 0.8]):
        source = runtime.replace("x = elementof(cartpow(reals, 4))", f"x = {literal(point)}")
        value = evaluate({"source": source, "bindings": ["result"]})[0]
        assert_allclose(value, runtime_oracle(point), rtol=3e-15, atol=0)
    example = (ROOT / "examples/spectrum-photometry.flatppl").read_text()
    source = example.replace("strengths = elementof(cartpow(nonnegreals, 2))", "strengths = [1.2, 0.7]")
    attenuated, source_counts, expected, score = evaluate({"source": source, "bindings": ["attenuated_flux", "source_counts", "expected_counts", "loglikelihood"]})
    template_values = 1e-12 * np.array([[1, 2], [2, 1], [1, 3]])
    expected_flux = (template_values @ [1.2, 0.7]) * np.exp(-np.array([0.2, 0.1, 0.05]) * np.log(10) / 2.5)
    expected_source = 100 * (np.array([[0.8, 0.5, 0], [0, 0.3, 0.7]]) @ (expected_flux * 1e-7 / (6.62607015e-34 * light / np.array([4e-7, 5e-7, 6e-7]))))
    assert_allclose(attenuated, expected_flux, rtol=2e-15, atol=0)
    assert_allclose(source_counts, expected_source, rtol=2e-15, atol=0)
    assert_allclose(expected, expected_source + [1, 2], rtol=2e-15, atol=0)
    assert_allclose(score, poisson.logpmf([70, 90], expected_source + [1, 2]).sum(), rtol=2e-14, atol=0)
    print(json.dumps({"checks": "redshift integrals and density Jacobian; photon/energy passbands; dust; unequal-area sky; runtime points; synthetic Poisson example", "example_counts": expected, "example_loglikelihood": score}))


if __name__ == "__main__":
    check()
