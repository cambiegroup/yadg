import pytest
import os
from tests.utils import (
    datagram_from_input,
    standard_datagram_test,
    compare_result_dicts,
)
import numpy as np
import json


def special_datagram_test(datagram, testspec):
    step = datagram["steps"][testspec["step"]]
    tstep = step["data"][testspec["point"]]


refvals = {
    "a": {"n": 4210.9762337, "s": 3.1576578, "u": " "},
    "b": {"n": 23349.8855948, "s": 16.2068347, "u": " "},
    "c": {"n": 71369.5831322, "s": 44.2360720, "u": " "},
    "d": {"n": 13119.4243886, "s": 7.23123702, "u": " "},
}


@pytest.mark.parametrize(
    "input, ts",
    [
        (
            {  # ts1 - dx file unzip, parse, method, integration from file
                "folders": ["."],
                "suffix": "txt",
                "parameters": {
                    "tracetype": "clarity.asc",
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
    print(ret)
    standard_datagram_test(ret, ts)
    special_datagram_test(ret, ts)


@pytest.mark.parametrize(
    "input",
    [
        (
            {  # ts0 - ch file parse, method, and integration
                "folders": ["."],
                "suffix": "CH",
                "parameters": {
                    "tracetype": "agilent.ch",
                    "calfile": "lc_calfile.json",
                },
            }
        ),
        (
            {  # ts1 - dx file unzip, parse, method, integration from file
                "folders": ["."],
                "suffix": "dx",
                "parameters": {
                    "tracetype": "agilent.dx",
                    "calfile": "lc_calfile.json",
                },
            }
        ),
    ],
)
def test_compare_raw_values(input, datadir):
    os.chdir(datadir)
    ret = datagram_from_input(input, "chromtrace", datadir)
    with open("yvals.json", "r") as infile:
        ref = json.load(infile)["traces"]
    for k, v in ret["steps"][0]["data"][0]["raw"]["traces"].items():
        for kk in ["t", "y"]:
            for kkk in ["n", "s"]:
                assert np.allclose(ref[k][kk][kkk], v[kk][kkk])


# def manual_test(datadir):
#     file = os.path.join(datadir, "dario.json")
#     with open(file, 'r') as f:
#         schema = json.load(f)
