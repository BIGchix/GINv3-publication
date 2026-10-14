"""Smoke: vendored SigPN engine loads without torch or xda/scBionet."""

from pathlib import Path


def test_scbionet_spn_modules_load():
    from gin_tool.subgraph_pipeline.spn_simulate import _load_scb_spn_modules

    mod_m, mod_s = _load_scb_spn_modules()
    assert hasattr(mod_m, "SPNModel")
    assert hasattr(mod_s, "SPNSimulator") or hasattr(mod_s, "simulate")


def test_node_list_and_seeds_present():
    root = Path(__file__).resolve().parents[1]
    assert (root / "nodeList" / "nodeList_TF.txt").is_file()
    assert (root / "seedLists/p53_pi3k_curated/genes.txt").is_file()
    assert (root / "seedLists/p53_pi3k_curated/active_tfs_rwrC_emp05.txt").is_file()
    assert (root / "data/GINv3.11a.sifi.gz").is_file()
    assert (root / "runtime/chip_gin_align/h1_wa01/allowlist_union.tsv").is_file()
