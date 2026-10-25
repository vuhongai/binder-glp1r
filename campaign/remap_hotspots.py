"""Round-2 hotspot remapping for the GLP1R campaign.

The round-2 input is a champion complex: chain A = binder (1..L), chain B =
target (renumbered L+1..L+N).  RFdiffusion's hotspot list must reference the
TARGET chain (B) with its renumbered positions.  Deriving the list from the
config's `reference_chain` + `hotspots` (as the TNF round-2 did) silently
points at the BINDER's own residues, and the partial diffusion runs unsteered.

This module maps by alignment, never a fixed offset: the champion's target
chain is aligned to the original target, each hotspot design-number is located
through the alignment, and the residue identity is verified at the mapped
position.  Any identity mismatch raises -- we never proceed with wrong
hotspots.

Self-test (CPU only, no GPU):

    python campaign/remap_hotspots.py --self-test data/targets/glp1r/receptor.pdb
"""
from __future__ import annotations

import os
import sys

ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}


def _three_to_one(name: str) -> str:
    name = name.strip().upper()
    if name in ONE:
        return ONE[name]
    # modified / unknown residue: fall back to the first letter
    return name[0] if name else "X"


def _read_chain(path: str, chain: str) -> list[tuple[int, str]]:
    """Ordered list of (residue_number, identity) for one chain."""
    out: list[tuple[int, str]] = []
    seen: set[int] = set()
    with open(path) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line[21] != chain:
                continue
            try:
                res = int(line[22:26])
            except ValueError:
                continue
            if res in seen:
                continue
            seen.add(res)
            out.append((res, _three_to_one(line[17:20])))
    out.sort(key=lambda p: p[0])
    return out


def _read_all_chains(path: str) -> dict[str, list[tuple[int, str]]]:
    chains: dict[str, list[tuple[int, str]]] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                chains.setdefault(line[21], None)
    for c in chains:
        chains[c] = _read_chain(path, c)
    return chains


def _nw_align(a: list[str], b: list[str]) -> list[tuple[int | None, int | None]]:
    """Global Needleman-Wunsch.  Returns (i, j) pairs: i indexes `a` (or None
    for a gap in a), j indexes `b` (or None for a gap in b)."""
    n, m = len(a), len(b)
    NEG = -10 ** 9
    S = [[NEG] * (m + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        S[i][0] = -i
    for j in range(m + 1):
        S[0][j] = -j
    for i in range(1, n + 1):
        ai = a[i - 1]
        row = S[i]
        prev = S[i - 1]
        for j in range(1, m + 1):
            s = prev[j - 1] + (1 if ai == b[j - 1] else -1)
            v = prev[j] - 1
            if v > s:
                s = v
            v = row[j - 1] - 1
            if v > s:
                s = v
            row[j] = s
    pairs: list[tuple[int | None, int | None]] = []
    i, j = n, m
    while i > 0 and j > 0:
        match = 1 if a[i - 1] == b[j - 1] else -1
        if S[i][j] == S[i - 1][j - 1] + match:
            pairs.append((i - 1, j - 1))
            i -= 1
            j -= 1
        elif S[i][j] == S[i - 1][j] - 1:
            pairs.append((i - 1, None))
            i -= 1
        else:
            pairs.append((None, j - 1))
            j -= 1
    while i > 0:
        pairs.append((i - 1, None))
        i -= 1
    while j > 0:
        pairs.append((None, j - 1))
        j -= 1
    pairs.reverse()
    return pairs


def compute_hotspots_contig(champion_pdb: str, target_pdb: str,
                            hotspots: list[int],
                            binder_chain: str = "A") -> tuple[str, list[tuple]]:
    """Map the original target's hotspot design-numbers onto the champion's
    target chain.  Returns (contig_string, details) where details is a list of
    (design_number, chain, champion_residue_number, identity).  Raises
    SystemExit on any identity mismatch or mapping failure."""
    if not os.path.isfile(champion_pdb):
        raise SystemExit(f"champion PDB not found: {champion_pdb}")
    if not os.path.isfile(target_pdb):
        raise SystemExit(f"target PDB not found: {target_pdb}")

    orig = _read_chain(target_pdb, "A")          # (design_number, identity)
    all_ch = _read_all_chains(champion_pdb)
    target_chains = [c for c in all_ch if c != binder_chain and all_ch[c]]
    if len(target_chains) != 1:
        raise SystemExit(
            f"expected exactly one target chain in {champion_pdb}, got "
            f"{target_chains} (binder chain = {binder_chain})")
    tchain = target_chains[0]
    champ = _read_chain(champion_pdb, tchain)    # (residue_number, identity)

    a = [id for _, id in orig]
    b = [id for _, id in champ]
    pairs = _nw_align(a, b)

    n_aligned = sum(1 for i, j in pairs if i is not None and j is not None)
    n_match = sum(1 for i, j in pairs
                  if i is not None and j is not None and a[i] == b[j])
    if n_match < n_aligned:
        raise SystemExit(
            f"target identity mismatch in {champion_pdb}: only {n_match}/"
            f"{n_aligned} aligned residues match the original target -- "
            "refusing to remap hotspots onto a different sequence")

    orig_index = {dn: idx for idx, (dn, _) in enumerate(orig)}
    imap: dict[int, int | None] = {}
    for i, j in pairs:
        if i is not None:
            imap[i] = j

    out: list[str] = []
    details: list[tuple] = []
    for h in hotspots:
        if h not in orig_index:
            raise SystemExit(
                f"hotspot {h} is not in the target (design 1-{len(orig)})")
        i = orig_index[h]
        j = imap.get(i)
        if j is None:
            raise SystemExit(
                f"hotspot {h} (design) has no aligned residue in the champion "
                "target chain")
        res_num = champ[j][0]
        if a[i] != b[j]:
            raise SystemExit(
                f"hotspot {h}: identity mismatch (original {a[i]} vs champion "
                f"{b[j]})")
        out.append(f"{tchain}{res_num}")
        details.append((h, tchain, res_num, a[i]))
    return "[" + ",".join(out) + "]", details


def _rewrite(line: str, chain: str, resseq: int) -> str:
    """Change the chain (col 21) and residue number (cols 22-25) of an ATOM
    line, preserving the file's (non-standard) column layout."""
    return line[:21] + chain + f"{resseq:4d}" + line[26:]


def _set_resname(line: str, resname: str) -> str:
    """Change the residue name (cols 17-19) of an ATOM line."""
    return line[:17] + resname.ljust(3) + line[20:]


def _self_test(target_pdb: str, hotspots: list[int], binder_len: int = 90) -> int:
    """Build a synthetic champion (binder_len-res binder + target renumbered
    +binder_len) by rewriting the real target lines, and confirm the remapper
    recovers the expected chain-B hotspots -- and that a wrong target identity
    is caught.  The expected contig is each hotspot design-number + binder_len."""
    import tempfile
    lines = [l for l in open(target_pdb) if l.startswith(("ATOM", "HETATM"))]
    if not lines:
        raise SystemExit(f"self-test: no ATOM lines in {target_pdb}")
    target = _read_chain(target_pdb, "A")
    if len(target) < binder_len + 1:
        raise SystemExit(f"self-test: target too short ({len(target)} res)")

    def write_champion(flip_resname_at: int | None) -> str:
        with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as fh:
            champ_path = fh.name
            # chain A: a binder_len-res binder (reuse the first target lines,
            # kept as chain A 1-binder_len -- the binder seq is irrelevant here)
            for l in lines[:binder_len]:
                fh.write(_rewrite(l, "A", int(l[22:26])))
            # chain B: the full target, renumbered +binder_len
            for l in lines:
                out = _rewrite(l, "B", int(l[22:26]) + binder_len)
                if flip_resname_at is not None and int(l[22:26]) == flip_resname_at:
                    out = _set_resname(out, "GLY")
                fh.write(out)
        return champ_path

    expected = "[" + ",".join(f"B{h + binder_len}" for h in hotspots) + "]"

    # happy path
    champ = write_champion(None)
    try:
        contig, details = compute_hotspots_contig(champ, target_pdb, hotspots)
        print(f"champion target chain : B (renumbered +{binder_len})")
        for h, c, r, ident in details:
            print(f"  design {h:3d}  ->  {c}{r}  ({ident})")
        print(f"hotspots_contig       : {contig}")
        print(f"expected              : {expected}")
        if contig != expected:
            print("SELF-TEST FAILED: contig mismatch")
            return 1
    finally:
        os.unlink(champ)

    # failure path: flip the identity of the first hotspot
    champ = write_champion(hotspots[0])
    try:
        try:
            compute_hotspots_contig(champ, target_pdb, hotspots)
            print("SELF-TEST FAILED: identity mismatch was not caught")
            return 1
        except SystemExit as exc:
            print(f"failure path OK       : raised ({str(exc).splitlines()[0]})")
        print("SELF-TEST PASSED")
        return 0
    finally:
        os.unlink(champ)


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] == "--self-test":
        if len(argv) < 3:
            raise SystemExit(
                "usage: remap_hotspots.py --self-test <target.pdb> "
                "<hotspot1,hotspot2,...> [binder_len]")
        hotspots = [int(x) for x in argv[2].split(",")]
        binder_len = int(argv[3]) if len(argv) > 3 else 90
        return _self_test(argv[1], hotspots, binder_len)
    raise SystemExit(__doc__)


if __name__ == "__main__":
    raise SystemExit(main())
