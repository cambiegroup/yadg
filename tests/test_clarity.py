import pytest
import os
from tests.utils import (
    datagram_from_input,
    standard_datagram_test,
    compare_result_dicts,
)


def special_datagram_test(datagram, testspec):
    step = datagram["steps"][testspec["step"]]
    tstep = step["data"][testspec["point"]]
    for k, v in testspec["peaks"].items():
        compare_result_dicts(tstep["derived"]["peaks"]["LC"][k]["A"], v)


refvals = {
    "a": {"n": 581.9376025004616, "s": 16.60366987726138, "u": " "},
    "b": {"n": 0.46264384587610796, "s": 0.05297980394927473, "u": " "},
    "c": {"n": 0.014848915735885246, "s": 0.0029636514936066827, "u": " "},
}


@pytest.mark.parametrize(
    "input, ts",
    [
        (
            {  # ts0 - ASCII file parse and integration from file
                "folders": ["."],
                "suffix": "txt",
                "parameters": {
                    "tracetype": "clarity.asc",
                    "calfile": "lc_calfile.json",
                },
            },
            {
                "nsteps": 1,
                "step": 0,
                "nrows": 1,
                "point": 0,
                "peaks": refvals,
            },
        ),
    ],
)
def test_datagram_from_lctrace(input, ts, datadir):
    os.chdir(datadir)
    ret = datagram_from_input(input, "chromtrace", datadir)
    standard_datagram_test(ret, ts)
    special_datagram_test(ret, ts)
