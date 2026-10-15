# GIN prior

| File | Role |
|------|------|
| `GINv3.11a.sifi.gz` | Frozen GINv3.11a edge list (default `--sifi`) |
| `GINv3.11a.sifi.gz.sha256` | Gzip sha256 + `uncompressed_sha256=…` |

The toolkit streams this file via transparent gzip I/O — **do not unpack** for
normal use. Fingerprints:

- gzip: `f9bb23c44c30d0dd6eab73225600085b378d1912e59f1ff0d8a3a0b37cfd0a58`
- uncompressed content: `192b2cb5a8147324fc41064a3f20c943a3d82c03972cff544e3d26219e7ace83`
  (matches lock JSON `gin.sha256`)

## Download

Shipped inside this repository under `data/` (`git clone` is enough).
GitHub Releases for this project **do not** re-attach the gzip as a binary
asset. Standalone download (pinned to tag `v0.2.0`):

- Browse: https://github.com/BIGchix/GINv3-publication/blob/v0.2.0/data/GINv3.11a.sifi.gz
- Raw: https://github.com/BIGchix/GINv3-publication/raw/v0.2.0/data/GINv3.11a.sifi.gz
- Sidecar: https://github.com/BIGchix/GINv3-publication/blob/v0.2.0/data/GINv3.11a.sifi.gz.sha256

Place the file at `data/GINv3.11a.sifi.gz` in a clone (or checkout) before
running `python scripts/reproduce_publication.py --steps prepare`.
