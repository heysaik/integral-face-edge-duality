# Exact checks for “Integral face-edge duality for knot triangulations”

These ancillary files accompany the preprint by Sai Kambampati. These finite checks supplement the proposed proofs; they do not establish universal statements or novelty.

## Reproduce the checks

Tested with Python 3.13. Install the three pinned packages in `requirements.txt`, then run:

```sh
python check_local_bridge.py
python check_affine.py
python verify.py
python check_rank_pairings.py
```

`check_local_bridge.py` verifies 48 local polynomial identities for both tetrahedron orientation signs. `check_affine.py` verifies the figure-eight affine constants, regular shapes, face isomorphism, and the need for a preferred-longitude qualification. `verify.py` reconstructs ten fixed ordered triangulations from `specimens.json` and checks:

- The full face-to-edge and preferred-longitude relation, using exact nullspaces.
- Nullity doubling and the nonsingular FAMED identity.
- Singular constraint-space and quadratic-form compatibility.
- Order-reversal Smith factors and the stronger transpose identity.
- Cyclic torsion of the face cokernel.
- Integral saturation of both lift matrices.
- The direct-sum decomposition of integer kernels, using saturated integer bases.
- Nullity relations modulo 2, 3, 5, 7, and 11, and odd-primary Smith factors.
- The index-two Neumann–Zagier row lattice.
- The determinant formula `abs(det B) = 2 * det(C)^2`.

The fixed specimens include face nullities 0, 1, and 2 in knot triangulations and integral homology solid tori with face determinants 3, 4, and 5. These latter manifolds are not asserted to be S³ knot exteriors. The face tables retain the actual ordering; an ordinary isomorphism signature alone would lose that information.

`verification-results.json` records the successful run. The longitude slope can vary with SnapPy's peripheral basis; it is recomputed by exact abelianization. The row uses the full signed-corner convention, with no division by two. Numerical shapes are never used to decide ranks, determinants, identities, or geometricity.

`check_rank_pairings.py` exhausts all 10,503 labelled ordered face-pairing tables through three tetrahedra. It filters for connected valid orientable one-torus ideal triangulations with first Betti number one, and checks that absence of edges of degree one or two implies a nonzero edge block. All 712 filtered tables, including 43 with no such low-degree edges and 31 strict-angle tables, pass. Results are recorded in `rank-pairing-results.json`. Counts include different labellings of the same manifold and are not counts of distinct manifolds.

## Exploratory evidence

The four older result files record separate, overlapping experiments:

| File | Scope |
|---|---|
| `probe-results.json` | 55 orders, rational ranks and nonsingular identities |
| `generalized-results.json` | 41 orders, singular constraints and forms |
| `bridge-results.json` | 29 orders, complete projected face relation |
| `move-results.json` | 54 strict-angle triangulations, 527 orders in a bounded Pachner-move determinant search |

These counts must not be added as if they were distinct specimens. Order enumeration identifies global reversal. The exploratory matrix code follows the published shape convention, including the exchange of SnapPy's prime and double-prime columns. Regina's strict-angle predicate is recorded separately and is not treated as a geometricity certificate.

## Mathematical sources

- Ben Aribi–Guilloux–Wong, [FAMED by computer](https://arxiv.org/abs/2512.17437), v1, definitions and Conjecture 4.1.
- Wong, [Generalized FAMED semi-geometric triangulations](https://arxiv.org/abs/2512.23198), v1, Definition 1.1.
- Wong, Section 2 of [Problems on Low-dimensional Topology, 2026](https://www.kurims.kyoto-u.ac.jp/~ildt/prob26.pdf), the kernel-correspondence and face-variable problems, with the figure-eight conventions used in `check_affine.py`.
- Garoufalidis–Hodgson–Hoffman–Rubinstein, [The 3D-index and normal surfaces](https://people.mpim-bonn.mpg.de/stavros/publications/normal.surfaces.3Dindex.pdf), integral normal-surface and lattice results used in the manuscript.

The scripts were written for this investigation. No upstream source files or third-party Python environments are included.
