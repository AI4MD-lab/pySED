"""Check installed pySED through its real CLI and independent analytic reference."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

import numpy as np


def check(output,entry_point="pysed",plot=True):
    output = Path(output).resolve()
    output.parent.mkdir(parents=True,exist_ok=True)
    environment = os.environ.copy()
    environment.pop("PYTHONPATH",None)  # Exercise the installed distribution.
    environment.update(MPLBACKEND="Agg",OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1")
    started = time.monotonic()

    def run(command,filename,cwd):
        result = subprocess.run(command,cwd=str(cwd),env=environment,
            stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=60)
        if output.is_dir():
            (output/filename).write_text(result.stdout,encoding="utf-8")
        if result.returncode:
            raise RuntimeError("Command failed: %r\n%s"%(command,result.stdout[-3000:]))
        return result.stdout

    run([sys.executable,str(Path(__file__).with_name("generate_demo.py")),"--output",str(output)],
        "generate.log",output.parent)
    run([sys.executable,"-m","pySED.prepare_vacancy_reference",
        "--perfect-basis","perfect_basis.in","--perfect-xyz","perfect.xyz",
        "--kept-rows","kept_rows.txt","--mapping-out","basis_vacancy.in",
        "--reference-out","perfect_reference.npz"],"prepare.log",output)
    if entry_point == "module":
        prefix = [sys.executable,"-m","pySED.main"]
    else:
        executable = Path(sys.executable).parent/(entry_point+(".exe" if os.name=="nt" else ""))
        if not executable.is_file():
            raise FileNotFoundError("Install pySED into this Python environment first: "+str(executable))
        prefix = [str(executable)]
    for input_file,log in [("input_SED.in","full.log"),("input_stream.in","stream.log")]:
        run(prefix+[input_file],log,output)
    with np.load(output/"expected_SED.npz",allow_pickle=False) as expected:
        scale = float(expected["SED"].max())
        if scale<=0:
            raise AssertionError("Reference cannot be identically zero")
        errors = {}
        for name in ["demo","stream"]:
            actual = np.loadtxt(output/(name+".SED"))
            error = float(np.max(np.abs(actual-expected["SED"]))/scale)
            if error>2e-12:
                raise AssertionError("SED differs from independent reference: %s"%error)
            np.testing.assert_allclose(np.loadtxt(output/(name+".Qpts")),expected["qpoints"],atol=1e-10,rtol=0)
            np.testing.assert_allclose(np.loadtxt(output/(name+".THz")),expected["frequency_THz"],atol=1e-7,rtol=0)
            errors[name] = error
    frequency = np.loadtxt(output/"demo.THz")
    spectrum = np.loadtxt(output/"demo.SED")
    peak = float(frequency[np.argmax(spectrum[:,1])])
    if abs(peak-7.8125)>1e-7:
        raise AssertionError("Injected peak was not recovered")
    if plot:
        run(prefix+["plot_SED.in"],"plot.log",output)
        if not (output/"demo-SED.png").is_file():
            raise AssertionError("Plotting CLI produced no PNG")
        settings = dict(line.split("=",1) for line in (output/"plot_SED.in").read_text().splitlines() if "=" in line)
        settings = {key.strip():value.strip() for key,value in settings.items()}
        peak_log = float(np.log(spectrum.max()))
        if not float(settings["colorbar_min"]) < peak_log <= float(settings["colorbar_max"]):
            raise AssertionError("Natural-log color limits hide the injected peak")
    origin = run([sys.executable,"-c","import pySED; print(pySED.__file__)"],"import.log",output).strip()
    report = dict(passed=True,entry_point=entry_point,imported_package=origin,
        normalized_max_errors=errors,injected_peak_THz=peak,plot_generated=plot,
        plot_peak_inside_color_range=True if plot else None,
        elapsed_seconds=time.monotonic()-started)
    (output/"validation.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,help="New output directory; defaults to a temporary directory")
    parser.add_argument("--entry-point",choices=["pysed","pySED","PYSED","module"],default="pysed")
    parser.add_argument("--no-plot",action="store_true")
    args = parser.parse_args()
    if args.output:
        result = check(args.output,args.entry_point,not args.no_plot)
    else:
        with tempfile.TemporaryDirectory(prefix="pysed-cli-") as temporary:
            result = check(Path(temporary)/"run",args.entry_point,not args.no_plot)
    print(json.dumps(result,indent=2))
