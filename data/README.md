# GIN prior

| File | Role |
|------|------|
| `GINv3.11a.sifi.gz` | Frozen GINv3.11a edge list (default `--sifi`) |
| `GINv3.11a.sifi.gz.sha256` | Gzip sha256 + `uncompressed_sha256=…` |

The toolkit streams this file via transparent gzip I/O — **do not unpack** for
normal use. Fingerprints:

- gzip: `f9bb23c4…`
- uncompressed content: `192b2cb5…` (matches lock JSON `gin.sha256`)
