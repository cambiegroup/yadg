"""
File parser for ClarityChrom ASCII export files (dat.txt).

This file format includes one timestep with a single trace in an ASCII file. It
contains a header section, and a sequence of Y datapoints.

Exposed metadata:
`````````````````

.. code-block:: yaml

    params:
      sampleid: !!str
      username: !!str
      datafile: !!str
      project: !!str

.. codeauthor:: Dario Cambie <dario.cambie@mpikg.mpg.de>
"""
import numpy as np

from yadg.dgutils.dateutils import str_to_uts


def process(fn: str, encoding: str, timezone: str) -> tuple[list, dict]:
    """
    ClarityChrom ASCII export file parser.

    One chromatogram per file with a single trace. A header section is followed by
    y-values for each trace.
    Method name is NOT available :(, nor is the detector name.
    They are assigned their numerical index in the file.

    Parameters
    ----------
    fn
        Filename to process.

    encoding
        Encoding used to open the file.

    timezone
        Timezone information. This should be ``"localtime"``.

    Returns
    -------
    ([chrom], metadata): tuple[list, dict]
        Standard timesteps & metadata tuple.
    """

    with open(fn, "r", encoding=encoding, errors="ignore") as infile:
        lines = infile.readlines()
    metadata = {"filetype": "claritychrom.datasc", "params": {}}
    chrom = {"fn": str(fn), "traces": {}}
    ytitle = None
    header_end = None
    samplerate = 0

    translate_key = {
        "Analyst": "username",
        "SampleID": "project",
        "Sample": "sampleid",
        "Data Original": "datafile",
    }

    for line_num, line in enumerate(lines):
        # Metadata
        for key in translate_key.keys():
            # note that the colon is needed as, e.g. "Sample Injected" starts with "Sample"
            if line.startswith(key+" :"):
                metadata["params"][translate_key[key]] = line.split(f"{key} :")[1].strip()

        # Chromatography parameters
        for key in ["Sample ID ", "Sample ", "Data Original "]:
            # note that the colon is needed as, e.g. "Sample Injected" starts with "Sample"
            if line.startswith(key+":"):
                k = key.lower().strip()
                chrom[k] = line.split(f"{key}:")[1].strip()

        if line.startswith("Sample Injected :"):
            chrom["uts"] = str_to_uts(
                line.split("Sample Injected :")[1].strip(),
                format="%d-%b-%y %I:%M:%S %p",
                timezone=timezone,
            )
        if line.startswith("Rate :"):
            assert (
                "per sec." in line
            ), f"clarityasc: Incorrect units for rate in file {fn}: {line}"
            samplerate_text = line.split(" : ")[1]
            samplerate = float(samplerate_text.split()[0])
        if line.startswith("Y Axis Title: "):
            ytitle = line[14:].strip()

        if ":" not in line:
            header_end = line_num
            break

    assert samplerate > 0, f"clarityasc: No samplerate found in file {fn}"
    assert isinstance(ytitle, str), f"clarityasc: No Y axis title in file {fn}"
    assert isinstance(header_end, int), f"clarityasc: Could not find any data in file {fn}"

    # X/Y units are not available in the header section, but are extracted from the data header.
    xunit, yunit = lines[header_end + 1].split("]	[")
    xunit = xunit.lstrip("[")
    yunit = yunit.strip().rstrip("]")

    points = len(lines) - (header_end + 2)
    x_axis, y_axis = zip(*[line.split() for line in lines[header_end + 2:]])
    x_axis = np.array(x_axis, dtype=float)

    # If time axis is in minute transform to second
    if xunit == "Min.":
        x_axis = x_axis * 60

    y_error = np.ones(points)
    # No info on Y error, assumed from units.
    if yunit == "mV":
        y_error = y_error * 0.01
    else:
        raise RuntimeError(f"Unknown y-axis error for {yunit}")

    y_axis = np.array(y_axis, dtype=float)
    x_error = np.ones(points) / samplerate

    assert len(x_axis) == len(y_axis) == len(x_error) == len(y_error)

    chrom["traces"][ytitle] = {
        "t": {
            "n": x_axis.tolist(),
            "s": x_error.tolist(),
            "u": "s",
        },
        "y": {
            "n": y_axis.tolist(),
            "s": y_error.tolist(),
            "u": yunit,
        },
        "id": 0,
        "data": [[x_axis, x_error], [y_axis, y_error]],
    }

    return [chrom], metadata
