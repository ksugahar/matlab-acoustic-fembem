# Readable Class Style

Gypsilab's source is valuable because the classes remain close to the
mathematics. This repository follows that spirit.

The goal is not to make the fastest MATLAB FEM/BEM solver. The goal is to make
a solver a student can open, read, and understand.

## Class Design Rules

- A class should represent one mathematical object:
  - mesh
  - finite-element space
  - boundary-element space
  - kernel
  - dense block
  - low-rank block
  - H-matrix tree
  - coupled FEM/BEM system
- Public properties should expose important mathematical data directly.
- Methods should be short enough to read during a lecture or code review.
- Constructors should only build the obvious object; expensive assembly should
  be a named method.
- Prefer explicit names such as `traceNodeIds`, `tetEdges`, `sourceWeights`,
  and `singleLayerCorrection`.
- Avoid clever vectorization when a simple loop teaches the formula better.
- Avoid hidden mutable caches unless they are clearly named and documented.
- Keep performance-oriented replacements behind separate functions or classes,
  not mixed into the teaching implementation.

## Good Shape

```matlab
space = H1Space(mesh);
K = space.stiffness();
M = space.mass();
```

The class file should make it obvious that `stiffness` assembles
`int grad(u) dot grad(v) dx`.

## Bad Shape

```matlab
space = OptimizedSpace(mesh, "CacheEverything", true);
K = space.op(17);
```

This may be fast, but it hides the mathematics and is not the design target.

## Performance Policy

NGSolve and NGSolve.BEM are the performance references. This repository can be
slower. If a faster implementation is needed later, keep it parallel to the
readable implementation and cross-validate it rather than replacing the
teaching code.

