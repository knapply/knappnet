# Country land borders (Natural Earth, 1:110m)

One-mode, undirected, spatial (node points and edge geometry), weighted
(border length): the 177 country units of Natural Earth's 1:110m "Admin 0
- Countries" layer and the pairs that share a land border, with the
border itself as edge geometry.

Natural Earth maps the world as it is administered de facto, so the
units include, for example, Kosovo, Somaliland, Northern Cyprus, Taiwan,
Western Sahara and Palestine, and treat Antarctica, Greenland, Puerto
Rico and the Falkland Islands as units; the layer's `type` column records
each unit's status. 21 units (islands and Antarctica) border nothing and
are isolates.

## Files

`nodes.csv` (177 rows), one per unit:

| column | meaning |
|---|---|
| `name` | Natural Earth's `ADMIN` name |
| `code` | Natural Earth's `ADM0_A3` code (unique) |
| `iso_a3` | ISO 3166-1 alpha-3 code, empty where Natural Earth gives none (France and Norway, whose polygons merge dependencies, and Kosovo, Northern Cyprus and Somaliland) |
| `continent`, `subregion` | Natural Earth's continent and UN subregion |
| `type` | `Sovereign country`, `Country`, `Dependency`, `Disputed`, `Indeterminate` or `Sovereignty` |
| `population`, `population_year` | Natural Earth's population estimate and its year |
| `longitude`, `latitude` | Natural Earth's label point for the unit (`LABEL_X`, `LABEL_Y`), WGS 84 (EPSG:4326) |

`edges.csv` (313 rows), one per pair of adjacent units:

| column | meaning |
|---|---|
| `from`, `to` | the two units' `code`s |
| `border_km` | length of the shared border, summed over great-circle arcs between consecutive vertices (mean Earth radius 6371.0088 km), rounded to 0.1 km; a strength-type weight |
| `geometry` | the shared border as a WKT `MULTILINESTRING` in WGS 84 longitude/latitude, one part per separate stretch (the Canada-United States border has two) |

Two units are adjacent when their polygons share at least one boundary
segment (identical consecutive vertex pairs, which the layer's
topologically consistent polygons guarantee for shared borders). Units
that touch only at a point are not adjacent. Lengths and shapes are those
of the generalized 1:110m geometry, so they are shorter and simpler than
the real borders (Canada-United States: 7660 km here).

## Source and build

Built by `data-raw/natural_earth.py` from
`geojson/ne_110m_admin_0_countries.geojson` in the Natural Earth vector
repository https://github.com/nvkelso/natural-earth-vector (commit
`ca96624a56bd078437bca8184e78163e5039ad19`, 2 June 2022, version
5.2.0-pre; SHA-256 recorded in the script). The same layer is published
at https://www.naturalearthdata.com/downloads/110m-cultural-vectors/.

## Licence

Public domain. Natural Earth's terms of use
(https://www.naturalearthdata.com/about/terms-of-use/): "All versions of
Natural Earth raster + vector map data found on this website are in the
public domain. You may use the maps in any manner, including modifying
the content and design, electronic dissemination, and offset printing."

## Citation

Natural Earth. Admin 0 - Countries, 1:110m cultural vectors, version
5.2.0-pre (2022). Made with Natural Earth, free vector and raster map data
at naturalearthdata.com.
