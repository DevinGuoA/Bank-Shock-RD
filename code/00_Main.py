"""Run the analyses without building the report or documentation."""

from pathlib import Path
import subprocess
import sys


### Paths ###

# Resolve all paths relative to the project, not the working directory.
ROOT = Path(__file__).resolve().parents[1]


### Run Analyses ###

for number in range(1, 5):
    script = ROOT / "code" / f"{number:02d}_Question_{number}.py"
    print(f"Running {script.name}", flush=True)
    # Reuse this interpreter and stop immediately if a question fails.
    subprocess.run([sys.executable, str(script)], cwd=ROOT, check=True)

print("All analyses completed. Compile the written report separately.")
