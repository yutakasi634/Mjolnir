+++
title = "g-BAOABLangevin"
weight = 2000
+++

# g-BAOABLangevin

`g-BAOABLangevin` integrator performs constant temperature simulation according to Langevin equation.

Unlike `BAOABLangevin`, it handles bond length constraints appropreately.

`g-BAOABLangevin` is developed in the following paper.

- Leimkuhler B, Matthews C. Proc. R. Soc. A. (2016)

## Example

```toml
[simulator]
integrator.type = "g-BAOABLangevin"
integrator.remove.translation = true
integrator.remove.rotation    = true
integrator.remove.rescale     = true
integrator.gammas = [
    {index = 0, gamma = 1.0},
    {index = 1, gamma = 1.0},
    # ...
]
```

## Input reference

Some of the other parameters, such as `delta_t`, are defined in [`[simulator]`]({{<relref "/docs/reference/simulators">}}) table.

- `type`: String
  - Name of the integrator. Here, it is `"g-BAOABLangevin"`.
- `remove`: Table (optional)
  - `translation` and `rotation`: Boolean
    - If `true`, it removes the total translation and rotation. Otherwise, it does nothing.
  - `rescale`: Boolean
    - If `true`, it rescales all the velocities to make kinetic energy constant.
  - By default, all the fields becomes `false`.
- `gammas`: Array of Tables
  - {{<katex>}}\gamma_i{{</katex>}} of the particles.

## Remarks

This feature is developed by contributor, [@yutakasi634](https://github.com/yutakasi634).
