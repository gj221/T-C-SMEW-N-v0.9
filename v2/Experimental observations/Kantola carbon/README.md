# Kantola harvest carbon observations

`maize_biomass_ERW_20.csv` is the unchanged block-level file from the
Illinois Data Bank dataset `10.13012/B2IDB-1917166_V1`. The dataset README
defines the three biomass columns as g C m^-2. Harvest biomass is grain for
the maize-soybean rotation.

`Load_Kantola_Harvest_Observations` calculates treatment-year means and
sample SEM directly from the five measured blocks. No dry-mass conversion
or assumed carbon fraction is used.
