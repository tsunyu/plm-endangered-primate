#!/usr/bin/env python3
"""Run a small SLiM nonWF cross-check for the inferred Ne schedule."""

from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from dfe_demography_common import (
    BASE,
    bin_frequencies,
    load_demography_mle,
    write_metadata,
)

OUT_ROOT = BASE / "output/method_validation/dfe_simulation"
SEED = 20260710


def build_slim_script(demography, n_loci: int, rep: int, out_file: Path) -> str:
    return f"""\
initialize() {{
    initializeSLiMModelType("nonWF");
    defineConstant("BURNIN", {demography.burn_in});
    defineConstant("TBOT", {demography.tbot_old});
    defineConstant("TREC", {demography.trecovery_old});
    defineConstant("TRECENT", {demography.trecent});
    defineConstant("NANC", {max(int(round(demography.nanc)), 50)});
    defineConstant("NBOT", {max(int(round(demography.nbot)), 50)});
    defineConstant("NREC", {max(int(round(demography.nrecover)), 50)});
    defineConstant("NCUR", {max(int(round(demography.ncur)), 50)});
    defineConstant("NLOC", {n_loci});
    defineConstant("ENDGEN", {demography.total_generations});
    defineConstant("OUTFILE", "{out_file.as_posix()}");
    initializeMutationRate(0.0);
    initializeRecombinationRate(0.0);
    initializeMutationType("m1", 0.25, "f", 0.02);
    m1.convertToSubstitution = F;
    initializeGenomicElementType("g1", m1, 1.0);
    initializeGenomicElement(g1, 0, NLOC - 1);
}}

1 early() {{
    sim.addSubpop("p1", NANC);
    for (pos in 0:(NLOC - 1)) {{
        sim.mutation(m1, pos, 0, 0, 0.01 + rgamma(1, 0.35, 0.01)[0]);
    }}
}}

1:ENDGEN late() {{
    gen = sim.generation;
    if (gen < BURNIN) {{
        p1.setSubpopulationSize(NANC);
    }} else {{
        since = gen - BURNIN;
        if (since < TBOT - TREC) {{
            p1.setSubpopulationSize(NBOT);
        }} else if (since < TBOT - TRECENT) {{
            p1.setSubpopulationSize(NREC);
        }} else {{
            p1.setSubpopulationSize(NCUR);
        }}
    }}
}}

ENDGEN late() {{
    freqs = sim.mutationFrequencies(p1, sim.mutations);
    if (size(freqs) == 0) freqs = c(0.0);
    writeFile(OUTFILE, freqs, append=F, sep="\\n");
    sim.simulationFinished();
}}
"""


def run_slim_crosscheck(replicates: int, n_loci: int, output: Path) -> pd.DataFrame:
    output.mkdir(parents=True, exist_ok=True)
    demography = load_demography_mle()
    rows = []
    for rep in range(replicates):
        freq_path = output / f"slim_rep_{rep}.freqs.txt"
        script = build_slim_script(demography, n_loci, rep, freq_path)
        with tempfile.NamedTemporaryFile("w", suffix=".slim", delete=False) as handle:
            handle.write(script)
            script_path = Path(handle.name)
        try:
            subprocess.run(["slim", str(script_path)], check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(exc.stderr or exc.stdout) from exc
        finally:
            script_path.unlink(missing_ok=True)
        if not freq_path.exists():
            raise FileNotFoundError(freq_path)
        freqs = np.loadtxt(freq_path, dtype=float)
        if freqs.ndim == 0:
            freqs = np.array([float(freqs)])
        bins = bin_frequencies(freqs)
        row = {"replicate": rep, "mean_daf": float(np.mean(freqs))}
        row.update({f"frac_{k}": v for k, v in bins.items()})
        rows.append(row)
    frame = pd.DataFrame(rows)
    frame.to_csv(output / "slim_crosscheck_summary.csv", index=False)
    write_metadata(
        output / "slim_crosscheck_metadata.json",
        {
            "replicates": replicates,
            "n_loci": n_loci,
            "seed": SEED,
            "slim_version": subprocess.check_output(["slim", "-version"], text=True).strip(),
            "note": "SLiM 5.1 nonWF cross-check using the same piecewise Ne schedule.",
        },
    )
    return frame


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUT_ROOT)
    parser.add_argument("--replicates", type=int, default=10)
    parser.add_argument("--loci", type=int, default=200)
    args = parser.parse_args()
    run_slim_crosscheck(args.replicates, args.loci, args.output)
    print(f"SLiM cross-check complete: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
