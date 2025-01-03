from argparse import ArgumentParser
from pathlib import Path
import pod5 as p5
from pod5.pod5_types import Calibration, Read
from pod5.tools.utils import collect_inputs
from tqdm import tqdm


def fix_calibration_scale(record: p5.ReadRecord) -> Read:
    """Fix the calibration scale for a read record by dividing by the digitisation"""
    if record.calibration.scale < 1:
        calibration = record.calibration
    else:
        calibration = Calibration(
            offset=record.calibration.offset,
            scale=record.calibration.scale / record.calibration_digitisation,
        )

    return Read(
        read_id=record.read_id,
        pore=record.pore,
        calibration=calibration,
        median_before=record.median_before,
        end_reason=record.end_reason,
        read_number=record.read_number,
        run_info=record.run_info,
        start_sample=record.start_sample,
        num_minknow_events=record.num_minknow_events,
        tracked_scaling=record.tracked_scaling,
        predicted_scaling=record.predicted_scaling,
        num_reads_since_mux_change=record.num_reads_since_mux_change,
        time_since_mux_change=record.time_since_mux_change,
        signal=record.signal,
    )


def argparser() -> ArgumentParser:
    parser = ArgumentParser(
        "fix_calibration_scale",
        description="Repairs the calibration scale value corrupted during file recovery.",
    )
    parser.add_argument(
        "-s",
        "--source",
        type=Path,
        default=Path.cwd(),
        help="Source directory of pod5 files. Does not search recursively.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path.cwd(),
        help="Output directory of pod5 files. Never overwrites files.",
    )
    parser.add_argument(
        "-x",
        "--suffix",
        type=str,
        default="",
        help="Suffix written to output pod5 file name as <source><.suffix>.pod5. Default: ``",
    )
    parser.add_argument(
        "-d",
        "--deep",
        action="store_true",
        help="Stops skipping files which appear normal.",
    )
    return parser


def main(source: Path, output: Path, suffix: str, deep: bool):
    """
    Fix calibration scale value in recovered pod5 files in source dir and write to
    output dir with optional filename suffix.
    """

    def out_path(src: Path, output: Path):
        if not suffix:
            return output / src.name
        return output / src.with_suffix(f".{suffix}.pod5").name

    print(f"Searching {source.resolve()} for pod5 files.")

    source_pod5s = collect_inputs([source], False, "*.pod5", threads=1)
    if not source_pod5s:
        raise ValueError(f"Found no pod5 files in: {source.resolve()}")
    else:
        print(f"Found {len(source_pod5s)} pod5 files to recover.")

    # Check output exists or make it
    if not output.exists():
        print(f"Creating output directory: {output.resolve()}")
        output.mkdir(parents=True)

    if not output.is_dir():
        raise NotADirectoryError(f"Output is not a directory: {output.resolve()}")

    # Check no overwrite
    existing = set(
        [out_path(s, output).name for s in source_pod5s if out_path(s, output).exists()]
    )
    if existing:
        raise FileExistsError(
            f"Found {len(existing)} output files which will not be overwritten."
        )

    # Iterate over all input files
    for s in tqdm(
        source_pod5s,
        desc="Processing Files",
        total=len(source_pod5s),
        dynamic_ncols=True,
        ascii=True,
        unit="file",
        position=0,
    ):
        with p5.Reader(s) as reader:
            # Skip empty files
            if reader.num_reads < 1:
                tqdm.write(f"zero reads in file: {s.resolve()}")
                continue

            if not deep:
                # Skip files that look normal
                first_record = next(reader.__iter__())
                if first_record.calibration.scale < 1:
                    tqdm.write(f"Calibration scale is normal in: {s.resolve()}")
                    continue

            with p5.Writer(out_path(s, output)) as writer:
                for record in tqdm(
                    reader,
                    "Fixing Reads",
                    total=reader.num_reads,
                    dynamic_ncols=True,
                    ascii=True,
                    unit="read",
                    smoothing=0,
                    mininterval=0.1,
                    position=1,
                    leave=False,
                ):

                    writer.add_read(fix_calibration_scale(record))

    print("Done")


if __name__ == "__main__":
    parser = argparser()
    args = parser.parse_args()
    main(**vars(args))
