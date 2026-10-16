"""
End-to-end self-calibration with real CASA: simulate a small VLA observation of a point source
with seeded phase errors and noise, then run phase self-calibration imaged with tclean.
"""
import os

import numpy as np
import pytest
from casatools import componentlist, measures, simulator, table

from snow.imaging import Tclean
from snow.selfcalibration import Phasecal

DIRECTION = ("J2000", "10h00m00s", "+30d00m00s")


def simulate_measurement_set(path: str) -> str:
    """Eight antennas, one spw, two hours of 60 s integrations of a 1 Jy point source."""
    rng = np.random.default_rng(1)
    east, north = rng.uniform(-2000, 2000, 8), rng.uniform(-2000, 2000, 8)
    me, sm, cl = measures(), simulator(), componentlist()
    cl.addcomponent(dir=" ".join(DIRECTION), flux=1.0, fluxunit="Jy", freq="5GHz", shape="point")
    cl.rename(path + ".cl")
    cl.close()
    sm.open(path)
    sm.setconfig(
        telescopename="VLA",
        x=east.tolist(),
        y=north.tolist(),
        z=[0.0] * 8,
        dishdiameter=[25.0] * 8,
        mount=["alt-az"] * 8,
        antname=["A{0}".format(i) for i in range(8)],
        padname=["P{0}".format(i) for i in range(8)],
        coordsystem="local",
        referencelocation=me.observatory("VLA")
    )
    sm.setspwindow(
        spwname="C",
        freq="5GHz",
        deltafreq="128MHz",
        freqresolution="128MHz",
        nchannels=1,
        stokes="RR LL"
    )
    sm.setfeed(mode="perfect R L")
    sm.setfield(sourcename="target", sourcedirection=me.direction(*DIRECTION))
    sm.setlimits(shadowlimit=0.001, elevationlimit="8deg")
    sm.setauto(autocorrwt=0.0)
    sm.settimes(
        integrationtime="60s",
        usehourangle=True,
        referencetime=me.epoch("utc", "2024/01/01/00:00:00")
    )
    sm.observe(sourcename="target", spwname="C", starttime="-1h", stoptime="1h")
    sm.predict(complist=path + ".cl")
    sm.close()
    me.done()
    corrupt(path, rng)
    return path


def corrupt(path: str, rng: np.random.Generator) -> None:
    """
    Antenna-based phase errors that self-calibration has to remove (a constant offset plus a slow
    drift), and thermal noise. Seeded numpy instead of the simulator's corruption, which is not
    reproducible from run to run.
    """
    tb = table()
    tb.open(path, nomodify=False)
    try:
        data = tb.getcol("DATA")  # (correlation, channel, row)
        antenna1, antenna2 = tb.getcol("ANTENNA1"), tb.getcol("ANTENNA2")
        hours = (tb.getcol("TIME") - tb.getcol("TIME").min()) / 3600.0
        offset = rng.uniform(-np.pi / 3, np.pi / 3, 8)
        drift, period = rng.uniform(-np.pi / 9, np.pi / 9, 8), rng.uniform(1.0, 3.0, 8)
        phase = offset[:, None] + drift[:, None] * np.sin(2 * np.pi * hours / period[:, None])
        phase[0] = 0.0  # the reference antenna
        rows = np.arange(len(hours))
        gains = np.exp(1j * (phase[antenna1, rows] - phase[antenna2, rows]))
        noise = 0.05 * (rng.standard_normal(data.shape) + 1j * rng.standard_normal(data.shape))
        corrupted = data * gains + noise
        # Imaging starts from CORRECTED_DATA, which must equal DATA before self-calibration
        for column in {"DATA", "CORRECTED_DATA"} & set(tb.colnames()):
            tb.putcol(column, corrupted)
    finally:
        tb.close()


def columns(path: str) -> list:
    tb = table()
    tb.open(path)
    try:
        return tb.colnames()
    finally:
        tb.close()


@pytest.fixture(scope="module")
def phasecal(tmp_path_factory):
    workdir = tmp_path_factory.mktemp("selfcal")
    # CASA writes scratch files to the working directory
    previous_cwd = os.getcwd()
    os.chdir(workdir)
    ms = simulate_measurement_set(str(workdir / "target.ms"))
    imager = Tclean(
        output=str(workdir / "target"),
        cell="2arcsec",
        M=128,
        N=128,
        niter=200,
        specmode="mfs",
        use_mask="user",  # no mask file: clean the whole image
        nsigma=3.0
    )
    phscal = Phasecal(
        visfile=ms,
        imager=imager,
        solint=["inf", "10min"],
        refant="A0",
        minblperant=3,
        want_plot=False
    )
    try:
        phscal.run()
    finally:
        os.chdir(previous_cwd)
    return phscal


def test_every_iteration_solves_and_images(phasecal):
    assert len(phasecal._caltables) == 2
    assert len(phasecal._psnr_history) == 3


def test_calibrated_data_and_model_are_written(phasecal):
    assert {"MODEL_DATA", "CORRECTED_DATA"} <= set(columns(phasecal.visfile))


def test_self_calibration_improves_the_image(phasecal):
    # The corruption is seeded, so this is deterministic: about 33 -> 55 -> 99
    psnr = phasecal._psnr_history
    assert psnr[0] < psnr[1] < psnr[2]
    assert psnr[2] > 2 * psnr[0]
