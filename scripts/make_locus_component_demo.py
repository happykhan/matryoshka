#!/usr/bin/env python3
"""Generate a compact visual acceptance artifact for the locus annotation vocabulary."""

from __future__ import annotations

import argparse
from pathlib import Path

from matryoshka.component_catalog import catalog_as_json, catalog_as_tsv
from matryoshka.detect import MGEFeature
from matryoshka.locus_map import to_locus_map_svg
from matryoshka.locus_table import to_locus_table_svg


def feature(
    element_type: str,
    name: str,
    start: int,
    end: int,
    strand: str = "+",
    family: str = "",
    **attributes: object,
) -> MGEFeature:
    return MGEFeature(
        element_type=element_type,
        family=family,
        name=name,
        start=start,
        end=end,
        strand=strand,
        attributes=attributes,
    )


def demo_hierarchy() -> list[MGEFeature]:
    region = feature("multiresistance_region", "multiresistance region", 100, 9800)
    integron = feature("integron", "class 1 integron", 450, 3900)
    cassette_array = feature("cassette_array", "cassette array", 1250, 2850)
    cassette_array.children = [
        feature("cassette", "aadA1a", 1260, 1900),
        feature("attC", "attC", 1850, 1900, strand="."),
        feature("cassette_remnant", "dfrA17", 1950, 2700, fragment=True),
        feature("attC", "attC#", 2660, 2700, strand="."),
    ]
    integron.children = [
        feature("IR", "IRi", 460, 484, strand="."),
        feature("integron_segment", "5'-CS", 500, 1200),
        feature("integrase", "intI1", 560, 980, strand="-"),
        feature("attI", "attI1", 1080, 1110, strand="."),
        feature("Pc_promoter", "Pc", 1140, 1170),
        cassette_array,
        feature("integron_segment", "3'-CS", 2851, 3700),
        feature("IR", "IRt", 3760, 3784, strand="."),
    ]
    tpu = feature("transposition_unit", "ISEcp1 TPU", 4100, 6750, family="ISEcp1_TPU")
    tpu.tsd_length = 5
    tpu.tsd_seq = "TATGA"
    tpu.children = [
        feature("IS", "ISEcp1", 4150, 5100),
        feature("captured_segment", "captured region", 5101, 6650),
        feature("AMR", "blaCTX-M-15", 5400, 6300),
        feature("IR", "IRalt", 6630, 6650, strand="."),
    ]
    region.children = [
        integron,
        tpu,
        feature("ISCR", "ISCR27#", 7000, 7700),
        feature("oriIS", "oriIS", 7010, 7030, strand="."),
        feature("terIS", "terIS", 7680, 7700, strand="."),
        feature("replicon", "IncFII", 7900, 8450),
        feature("oriV", "oriV", 8500, 8530, strand="."),
        feature("ncRNA", "RNAI", 8650, 8950, strand="-"),
        feature("direct_repeat", "DR", 9100, 9104, strand=".", sequence="AATGC"),
        feature("unknown_fragment", "unresolved", 9300, 9600),
    ]
    return [region]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    roots = demo_hierarchy()
    (args.out / "locus-component-vocabulary.svg").write_text(
        to_locus_map_svg(roots, 10_000, "locus component vocabulary"),
        encoding="utf-8",
    )
    (args.out / "locus-component-table.svg").write_text(
        to_locus_table_svg(roots, "locus component vocabulary"),
        encoding="utf-8",
    )
    (args.out / "locus-component-catalog.tsv").write_text(
        catalog_as_tsv(),
        encoding="utf-8",
    )
    (args.out / "locus-component-catalog.json").write_text(
        catalog_as_json(),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
