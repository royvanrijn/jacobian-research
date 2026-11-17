# Lean formalizations

The developments in this directory are separate Lake packages with independent
targets and manifests. They use the same Lean and Mathlib release, so their
remote dependencies are stored once in the ignored directory
`formal/.lake/packages`. Each `lakefile.toml` selects that directory through
the relative setting:

```toml
packagesDir = "../.lake/packages"
```

Existing commands should still be run from the individual package directory.
For example:

```bash
cd formal/gvc
lake exe cache get
lake build
```

Running `lake update` in a package materializes missing dependencies in the
shared directory while retaining that package's independent lock manifest.
Avoid concurrent `lake update`, `lake exe cache get`, or `lake clean` commands
across these packages because those commands mutate the shared dependency
tree. If a package moves to a different Lean or Mathlib pin, give that package
a version-specific `packagesDir` instead of sharing incompatible dependencies.
