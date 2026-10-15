# Examples

| Path | Purpose |
|------|---------|
| `01_custom_context.sh` | Toolkit: your seeds → context subgraph (uses `.sifi.gz` or plain) |
| `02_publication_stage_a.sh` | Publication: **recompute** Stage A, print sha vs lock |
| `03_sigpn_demo.sh` | Toy subgraph → SigPN export + short MC (**no** GIN prior) |
| `run_sigpn_demo.py` | Python entry used by `03_sigpn_demo.sh` |
| `demo_sigpn/` | Bundled toy edges + place levels for the SigPN demo |
| `allowlist_stub/` | Tiny ChIP allowlist **format** demo (not real peaks) |

Docs: `docs/USER_GUIDE.md` · `docs/CHIP_ALLOWLIST.md` · `docs/SIGPN_CUSTOM_SIMULATIONS.md`

```bash
python scripts/reproduce_publication.py --steps prepare
bash examples/01_custom_context.sh
bash examples/03_sigpn_demo.sh
N_JOBS=20 bash examples/02_publication_stage_a.sh
```
