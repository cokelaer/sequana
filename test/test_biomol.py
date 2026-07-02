import pytest

from sequana.biomol import *


def test_biomol():

    assert compute_melting_temperature_wallace_rule("AGTGGGTTCCGT") == 38
    assert pytest.approx(compute_melting_temperature_salt_adjusted("AGTGGGTTCCGT")) == 39.05

    # lower-case (soft-masked) input must give the same result
    assert compute_melting_temperature_wallace_rule("agtgggttccgt") == 38
    assert pytest.approx(compute_melting_temperature_salt_adjusted("agtgggttccgt")) == 39.05
