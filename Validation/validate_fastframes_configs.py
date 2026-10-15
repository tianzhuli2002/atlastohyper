#!/usr/bin/env python3
"""Parse FastFrames configs without launching event processing."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


FASTFRAMES_ROOT = Path("/home/leet/FF/fastframes")


def validate(paths: list[Path]) -> list[dict]:
    if not FASTFRAMES_ROOT.is_dir():
        raise FileNotFoundError(f"FastFrames source is not available: {FASTFRAMES_ROOT}")

    sys.path[:0] = [
        str(FASTFRAMES_ROOT),
        str(FASTFRAMES_ROOT / "python"),
        str(FASTFRAMES_ROOT / "python" / "ConfigReaderModules"),
    ]
    from BlockReaderCommon import set_paths

    set_paths()
    from CommandLineOptions import CommandLineOptions
    from ConfigReader import ConfigReader

    results = []
    for config_path in paths:
        config_path = config_path.resolve()
        if not config_path.is_file():
            raise FileNotFoundError(f"FastFrames config is missing: {config_path}")

        # ConfigReader consults FastFrames' CLI singleton while resolving sample
        # settings. Reset it so each requested config receives its own path.
        CommandLineOptions._instances.clear()
        sys.argv = ["validate_fastframes_configs.py", "--config", str(config_path), "--step", "n"]
        reader = ConfigReader(str(config_path))
        filelist = Path(reader.block_general._input_filelist_path)
        sumweights = Path(reader.block_general._input_sumweights_path)
        if not filelist.is_file():
            raise FileNotFoundError(f"{config_path}: filelist does not exist: {filelist}")
        if not sumweights.is_file():
            raise FileNotFoundError(f"{config_path}: sum-of-weights file does not exist: {sumweights}")
        results.append(
            {
                "config": str(config_path),
                "parse": "passed",
                "regions": list(reader.regions),
                "samples": list(reader.samples),
                "input_filelist_path": str(filelist),
                "input_filelist_exists": True,
                "input_sumweights_path": str(sumweights),
                "input_sumweights_exists": True,
            }
        )
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("configs", nargs="+", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    result = {"configs": validate(args.configs)}
    output = json.dumps(result, indent=2) + "\n"
    print(output, end="")
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(output)


if __name__ == "__main__":
    main()
