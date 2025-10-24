# Convertable

Conversion trait for:
- `alloy-primitives`
- `malachite`
- `bigdecimal`


### Features:
```toml
[features]
default = ["full"]

full = ["malachite", "alloy", "bigdecimal"]
malachite = ["dep:malachite"]
alloy = ["dep:alloy-primitives"]
bigdecimal = ["dep:bigdecimal"]
```