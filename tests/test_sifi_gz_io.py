"""Smoke: plain .sifi and .sifi.gz yield identical edge streams."""

from __future__ import annotations

import gzip
import tempfile
from pathlib import Path

from gin_tool.io import collect_unique_node_names, iter_sifi_edges, open_sifi_text


_TSV = (
    "TP53\tMDM2\tpositive\n"
    "MDM2\tTP53\tnegative\n"
    "ATM\tTP53\tpositive\n"
)


def test_open_sifi_text_plain_and_gz_match() -> None:
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        plain = root / "toy.sifi"
        gz = root / "toy.sifi.gz"
        plain.write_text(_TSV, encoding="utf-8")
        with gzip.open(gz, "wt", encoding="utf-8") as fh:
            fh.write(_TSV)

        with open_sifi_text(plain) as a, open_sifi_text(gz) as b:
            assert list(a) == list(b)

        assert list(iter_sifi_edges(str(plain))) == list(iter_sifi_edges(str(gz)))
        assert collect_unique_node_names(str(plain)) == collect_unique_node_names(
            str(gz)
        )
        assert collect_unique_node_names(str(gz)) == {"TP53", "MDM2", "ATM"}
