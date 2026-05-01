# Finding 1: Sample-Selection Mismatch Between the Notebook and the Standalone 90CL Script

## Summary

The current standalone script does not select posterior samples in exactly the same way as the original XMM notebook.

This is not a general disagreement about the pressure-profile model itself. The issue is narrower:

- the notebook and the standalone script are constructing the profile ensemble from slightly different subsets of the posterior samples
- because the percentile bands are computed from that ensemble, even a one-sample difference can shift the plotted 5th, 50th, and 95th percentile curves

This is the main concrete logic divergence found between the original notebook and the current combined script.

## Where the Difference Comes From

### Notebook behavior

In the notebook, the posterior samples are first obtained from:

- `samples, weights, logl, _ = sampler.posterior()`

Later, the profile ensembles are built with a loop over 1000 samples, and each iteration accesses:

- `samples[-i]`

The important detail is that Python treats `-0` as `0`.

That means the notebook is not selecting the last 1000 samples as a contiguous block.

Instead, it selects:

- `samples[0]`
- `samples[-1]`
- `samples[-2]`
- `samples[-3]`
- ...
- `samples[-999]`

So the notebook sample set is:

```python
samples[0], samples[-1], samples[-2], ..., samples[-999]
```

### Standalone script behavior

The standalone script loads the flattened sample array from the saved state and returns:

```python
samples[-sample_count:]
```

For `sample_count = 1000`, this selects:

```python
samples[-1000], samples[-999], ..., samples[-1]
```

So the standalone script sample set is:

- `samples[-1000]`
- `samples[-999]`
- ...
- `samples[-1]`

## Exact Difference in Indexing

The two sample-selection rules differ by one sample:

- the notebook includes `samples[0]`
- the standalone script includes `samples[-1000]`

Everything else in the chosen 1000-sample set overlaps.

That may look minor, but it is still a real difference in logic.

## Why This Affects the Plots

The XMM uncertainty bands are not drawn from a closed-form error formula. They are built empirically from an ensemble of profiles generated from the chosen posterior samples.

The workflow is:

1. choose a set of posterior samples
2. convert each sample into a physical pressure profile
3. evaluate those profiles on a radius grid
4. take the 5th, 50th, and 95th percentiles at each radius

That means the plotted median line and shaded band depend directly on which sample set is used.

Even changing one sample can move the percentile boundaries, especially the lower and upper tails:

- the 50th percentile can move slightly
- the 5th percentile can move more noticeably if the replaced sample is relatively extreme
- the 95th percentile can also shift if the swapped sample sits near the upper envelope

For the current saved state, the difference is small for most outputs, but it is not zero. The largest measured effect was in the west-subcluster lower percentile curve, where the relative change was about 7.7%.

So the key point is:

- this is not just a cosmetic indexing detail
- it changes the ensemble used to define the plotted uncertainty region
- therefore it can change the plotted XMM band

## Why This Likely Happened

The most likely explanation is that the notebook author intended one of these two behaviors:

### Possibility 1: use the last 1000 samples

If that was the intention, the correct code would have been:

```python
samples[-1000:]
```

### Possibility 2: step backward through the last 1000 samples one by one

If that was the intention, the loop should have used:

```python
samples[-(i + 1)]
```

instead of:

```python
samples[-i]
```

Because `i = 0` in the first iteration, `samples[-i]` becomes `samples[0]`, which is almost certainly not the intended “last sample” behavior.

So the most plausible interpretation is that the notebook contains a small indexing quirk rather than a deliberate statistical choice.

## Why This Matters for Reproducibility

There are two legitimate goals, but they are different:

### Goal A: exact notebook reproduction

If the goal is to reproduce the notebook literally, then the standalone script should intentionally match the notebook's sample-selection rule, including the `samples[0]` inclusion.

### Goal B: cleaned-up standalone implementation

If the goal is to implement the likely intended logic, then using the last 1000 samples as a contiguous block is more natural and easier to explain.

The current issue is that the standalone script is close to the notebook, but not exactly identical, while the comments imply that it is mirroring the notebook behavior.

That makes the current implementation ambiguous from a reproducibility standpoint.

## Practical Interpretation

At the moment, there is no evidence that this difference creates a catastrophic error in the plotted profiles.

What it does mean is:

- the combined 90CL script is not an exact drop-in reproduction of the notebook profile-construction logic
- any agreement between the notebook-style output and the standalone output is approximate rather than exact
- future readers of the code could easily assume the two are identical when they are not

## Recommended Resolution

The project should make an explicit choice between two options.

### Option 1: preserve notebook fidelity

If exact reproduction matters most, then the standalone script should deliberately match the notebook indexing rule and document that choice.

### Option 2: preserve the likely intended logic

If clarity and maintainability matter most, then the standalone script should continue using the last 1000 samples as a contiguous slice, but the comments should say clearly that this is an intentional cleanup rather than a literal transcription.

## Recommended Follow-Up

The cleanest follow-up would be:

1. add a small helper dedicated to sample selection
2. give that helper a name that makes the choice explicit
3. add a regression test that locks in whichever interpretation the project decides to keep

That would remove the ambiguity and make future comparisons between the notebook and the standalone script easier to interpret.