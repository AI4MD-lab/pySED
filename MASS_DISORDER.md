# Opt-in substitutional mass disorder

allow_mass_disorder = 1 requires an explicit complete host reference.
Mapping rows retain the physical host basis index but record the actual atom mass.
For each host basis b and Cartesian direction, compute the coherent amplitude
sum_i sqrt(m_i) v_i(t) exp(i q.R_parent_cell_i), then its temporal Fourier power.
The implementation multiplies v_i by sqrt(m_i/m_reference_b) before the sum and
multiplies the power by m_reference_b. The reference mass cancels analytically.
Vacant sites contribute zero and normalization uses the complete host cell count.
This includes Nb-Ti cross terms; grouping them into separate species powers would
drop those terms. Partial output is grouped by distinct reference basis masses and
Cartesian directions, as in pySED's existing output convention. For NbFeSb, three
distinct reference masses give nine channels; the Nb channel also includes Ti.
Different parent bases with identical reference masses share a mass-type channel.
These outputs are not pure-species PDOS.

Default 0 retains the vacancy-only mass consistency checks and original spectra.
Uniform masses reproduce the old formula. Streaming and full-read paths both weight
the actual per-row mass. This is an eigenvector-free mass-weighted host-q projection,
not an eigenvector-resolved unfolding or proof of defect-potential accuracy.

Mass-weighted normal-mode coordinates are described in the primary thesis:
https://jasonlarkin.github.io/cmuthesis_jlarkin.pdf (section 2.2).
The extension to individual site masses here defines the disordered projection.

For a real substitution, update only the actual mass in column four of the
surviving-atom mapping; retain its original cell and parent-basis indices. For
example, Ti on a Nb site keeps the Nb basis index and records 47.867 amu instead
of 92.90637 amu. Do not change the mass of that basis in the complete parent NPZ.
The default remains `allow_mass_disorder = 0`.
