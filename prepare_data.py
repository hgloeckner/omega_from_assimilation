# %%
import xarray as xr
import pydropsonde.circles as circles
import numpy as np
from orcestra import get_flight_segments

meta = get_flight_segments()
# %%
original_vars = [
    f"{var}_{kind}"
    for var in ["Uwind", "Vwind", "Hum", "Temp", "P"]
    for kind in ["Obs", "AN", "FG"]
] + ["launch_datetime"]
new_vars = [
    f"{var}_{kind}"
    for var in ["u", "v", "q", "ta", "p"]
    for kind in ["obs", "an", "fg"]
] + ["launch_time"]
# %%
ds = xr.open_dataset("DS_IFS_assimilatedODB_All.nc")
ds = ds.where(ds.launch_datetime != "--", drop=True).rename(
    {orig: new_var for orig, new_var in zip(original_vars, new_vars)}
)
lev4 = xr.open_dataset(
    "ipfs://QmP73Kosem4exJcZXxG8vpN4YLqaepoZSPWwnQ9N1xffus/Level_4/PERCUSION_Level_4.zarr",
    engine="zarr",
)
lev4 = (
    lev4.assign_coords(
        level=("altitude", np.arange(lev4.sizes["altitude"])),
    )
    .swap_dims({"altitude": "level"})
    .interp(level=ds.level)
    .swap_dims({"level": "altitude"})
)


ds = ds.assign(launch_time=ds.launch_time.astype(np.datetime64))
# %%


# %%


def get_circle_id_from_sonde_id(lev4, sonde_id):
    sonde_idx = (
        lev4.swap_dims({"sonde": "sonde_id"}).get_index("sonde_id").get_loc(sonde_id)
    )
    circle_end = lev4.sondes_per_circle.cumsum()
    cidx = np.argmin(np.abs(sonde_idx - circle_end.values))
    if circle_end[cidx] > sonde_idx:
        return cidx
    else:
        return cidx + 1


def get_sonde_ids_from_l4(ds, lev4):
    ids = []
    wrong_sondes = []
    for sonde in ds.profile.values:
        lat = ds.sel(profile=sonde).latitude.values
        lon = ds.sel(profile=sonde).longitude.values
        time = ds.sel(profile=sonde).launch_time.values
        proxy_l4 = lev4.swap_dims({"sonde": "launch_time"}).sel(
            launch_time=slice(
                time - np.timedelta64(5, "m"), time + np.timedelta64(5, "m")
            )
        )
        if proxy_l4.sizes["launch_time"] == 0:
            wrong_sondes.append(sonde)
        else:
            idxmin = np.argmin(
                np.abs(lat - proxy_l4.aircraft_latitude.values)
                + np.abs(lon - proxy_l4.aircraft_longitude.values)
            )
            ids.append(
                ds.sel(profile=sonde)
                .assign(
                    sonde_id=proxy_l4.sonde_id.values[idxmin],
                    u_orc=("level", proxy_l4.u.values[idxmin]),
                    v_orc=("level", proxy_l4.v.values[idxmin]),
                    q_orc=("level", proxy_l4.q.values[idxmin]),
                    ta_orc=("level", proxy_l4.ta.values[idxmin]),
                    p_orc=("level", proxy_l4.p.values[idxmin]),
                    circle_id=lev4.isel(
                        circle=get_circle_id_from_sonde_id(
                            lev4, proxy_l4.sonde_id.values[idxmin]
                        )
                    ).circle_id.values,
                )
                .set_coords(
                    ["sonde_id", "circle_id", "longitude", "latitude", "launch_time"]
                )
            )

    return ids, wrong_sondes


# %%
ids, wrong_sondes = get_sonde_ids_from_l4(ds, lev4)
# %%
circle_ds = xr.concat(ids, dim="profile")
circle_ds.to_netcdf("ds_assimilated_and_orcestra.nc")
# %%

import matplotlib.pyplot as plt
import seaborn as sns
import cartopy.crs as ccrs
import cartopy.feature as cf

fig, ax = plt.subplots(
    subplot_kw={"projection": ccrs.PlateCarree()},
)
ax.set_extent([-70, -10, -2, 18], crs=ccrs.PlateCarree())
ax.add_feature(cf.LAND, facecolor="lightgray")
ax.add_feature(cf.COASTLINE, edgecolor="black")
for sonde in wrong_sondes:
    ax.plot(
        ds.sel(profile=sonde).longitude.values,
        ds.sel(profile=sonde).latitude.values,
        marker="o",
        markersize=2,
        color="red",
    )
for sonde in ids:
    ax.plot(
        sonde.longitude.values,
        sonde.latitude.values,
        marker="o",
        markersize=2,
        color="blue",
    )
