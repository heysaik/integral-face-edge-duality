# Integral face-edge duality for knot triangulations

**Sai Kambampati**  
[sai@kambampati.studio](mailto:sai@kambampati.studio)

[Read the manuscript (PDF)](face-cochains-famed.pdf) · [LaTeX source](face-cochains-famed.tex) · [Computational checks](anc/README.md)

## Status

This mathematical preprint contains proposed proofs and is circulated for critical scrutiny. Partial informal specialist feedback has been received. Full independent verification of the arguments and their novelty, and formal peer review, remain incomplete. The current revision is dated October 2, 2026.

The release records the exact manuscript and supporting files made public at that time. It does not establish mathematical correctness, priority over other work, or acceptance by arXiv or a journal. The preprint was submitted to arXiv in math.GT on October 2, 2026 and is awaiting processing and moderation. No public arXiv identifier or DOI has been assigned yet.

## Revision history

- **v2, October 2, 2026:** clarifies definitions, matrix and peripheral conventions, and the translation from cited foundational results to the manuscript's coordinates. The theorem claims and computational scripts are unchanged.
- **[v1, September 18, 2026](https://github.com/heysaik/integral-face-edge-duality/releases/tag/v1):** initial public release of the manuscript dated September 16, 2026. This snapshot remains available.

## Scope

The paper proposes a face-cochain correspondence for ordered ideal triangulations of knot exteriors in the three-sphere, using the preferred longitude. It relates the face matrices to the gluing matrices, including singular cases, and develops integral kernel and cokernel consequences. Under a strict angle structure, it gives a proposed proof of Wong's generalized FAMED conjecture.

The determinant-one assertion, the extra longitude-meridian conditions, and the full Andersen-Kashaev volume conjecture remain open in this work. Exact computations check finite examples and conventions; the general claims depend on the mathematical arguments in the manuscript.

## Files and reproduction

- `face-cochains-famed.pdf`: the 26-page manuscript, including references and the final appendix on AI use and verification status.
- `face-cochains-famed.tex`: its self-contained LaTeX source; the figure is drawn directly in LaTeX.
- `anc/`: exact-arithmetic scripts, pinned dependencies, input specimens, and recorded results.
- `SHA256SUMS.txt`: file checksums for this snapshot.
- `CITATION.cff`: citation metadata.

For the computational checks, follow [anc/README.md](anc/README.md). They were tested with Python 3.13 and the package versions listed in `anc/requirements.txt`.

The PDF was built with Tectonic. To build it locally, run:

```sh
tectonic face-cochains-famed.tex
```

To verify the downloaded snapshot on macOS, run:

```sh
shasum -a 256 -c SHA256SUMS.txt
```

A local LaTeX rebuild may differ byte-for-byte from the archived PDF because of the TeX environment and generated metadata. The checksum identifies the published PDF itself.

## Feedback

Focused corrections, counterexamples, and references to prior work are welcome through GitHub issues or the author's email. Please identify the relevant statement and hypothesis when reporting a mathematical issue. Suggested starting points are the preferred-longitude argument in Section 5, the generalized-FAMED conclusion in Section 8, and the integral statements in Sections 10-12.

arXiv endorsement is separate from mathematical review.

## Citation

Sai Kambampati, *Integral face-edge duality for knot triangulations*, preprint, GitHub release v2 (October 2, 2026).

Please cite a specific release or commit when referring to this version. Later corrections will be recorded in new commits and releases.
