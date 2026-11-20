# fig5a_TIL_zones_map.ipynb (new)

Origin: `_plot_til_zones_ax` / `_til_percent_by_sample` of EFS Analysis_Repo/Notebooks/TIL_Niche_Comparisons_v2_clean.ipynb cell 26 (published style, Spencer Watson); `plot_visium_til` of pipelines/NicheNet_v3.ipynb cells 184–185.

Changes:
- Sections set explicitly: ER_G_631a_vis (TIL-high) and CR_G_310a_vis (TIL-low), identified by matching the published TIL+ islands; the facet cell's automatic top/bottom selection picks ER_G_740a/736a and CR_G_322a/310a with the frozen zones (printed).
- `_plot_til_zones_ax` takes `adata` as an argument; coordinates from obs only; returns coords for the new `add_scale_bar` (1 mm, pixel size from the 100 µm spot pitch).
- `_til_percent_by_sample`: spot count aggregates `("_til", "size")` (the original `(key, "count")` raises KeyError with as_index=False in pandas 2.2).
- `plot_visium_til` verbatim except: no `plt.show()` inside (so the figure can be saved), NearestNeighbors import at top.
- Data: `DataCenter.load_spatial(list_sample_id=..., load_gex=False)` per centre (AllDataCenters.load_spatial with a list fails for centres without listed sections) + frozen `config.TIL_ZONES_CSV`.

Package helpers: `config` (paths), `output`.
Verification: TIL+ counts 660 (ER_G_631a) / 7 (CR_G_310a); maps match the published Fig 5a (published panels are crops).
