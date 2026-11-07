# Bounded-memory projection of long trajectories

The optional streaming projection keeps the coherent time signals for all q-points
and reads atomic velocities in short I/O blocks:

```text
sed_projection_block_size = 128
compress_block_size = 128
num_splits = 1
```

Both block-size controls default to zero. Zero projection block size selects the
existing full-trajectory calculation; zero compression block size retains upstream
automatic block selection. A positive projection block size selects the streamed
q-batched path. It does not use q-point multiprocessing; BLAS threading is controlled
through the numerical-library environment rather than `max_cores`.

For every basis and direction, each I/O chunk fills a slice of the complete
coherent amplitude time series by a matrix product:

```python
phase = exp(1j * reference_origins @ qpoints.T)
signal[start:stop, :, b, alpha] = velocity_chunk[:, ids, alpha] @ phase[cell_ids]
```

Only after all chunks for the split have been assembled is an FFT performed along
time. Squared FFT amplitudes, basis masses, partial-spectrum grouping and upstream
normalization are then applied normally. The powers of separate I/O chunks are
never averaged. Consequently, I/O chunking does not shorten the FFT window or
change its frequency resolution. `num_splits` continues to control the physical
FFT-window length and averaging independently.

For `F` frames, `A` atoms, `Q` q-points, `B` basis sites and chunk size `K`, atomic
trajectory working memory scales with `K*A*3` rather than `F*A*3`. Complex signals
still scale with `F*Q*B*3`; streaming is not an unlimited-memory guarantee.
Partial SED arrays also scale with frame count and q count.

For 20,000 frames, 24,000 parent atoms and 33 q-points, a full float64 velocity
array alone is 11.52 GB. A 128-frame velocity chunk is 73.728 MB and the complex
signal array is 95.04 MB for three basis sites. These are array-size estimates,
not measured total process RSS; see the benchmark resource logs for measurements.

The GPUMD text reader uses NumPy's numeric parser for each frame, extracting the
six position/velocity columns. It checks atom count, supported property order,
shape and finite values. The NumPy minimum is raised to 1.16 for
[`max_rows`](https://numpy.org/doc/1.16/reference/generated/numpy.loadtxt.html).
The existing LAMMPS input layout is retained.

New HDF5 caches are marked incomplete until compression finishes. An explicitly
incomplete cache is rejected during lattice initialization rather than silently
reusing its preallocated but unfilled rows. Older caches without this attribute
remain compatible.

Tests compare streamed and full calculations with I/O chunk lengths of 1, 7 and
64 frames, uneven chunk endings, two physical trajectory splits, missing
sublattices and the legacy complete-crystal phase reference. Additional tests
verify consecutive GPUMD frame parsing and reject an unsupported column order.
