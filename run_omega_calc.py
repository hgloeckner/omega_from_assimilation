# %%
import xarray as xr
import pydropsonde.circles as circles
import numpy as np
from orcestra import get_flight_segments
from attrs import attr_dict, circle_attr_dict

# %%
meta = get_flight_segments()
ds = xr.open_dataset("ds_assimilated_and_orcestra.nc").load()
details = [
    segment
    for flight in meta["HALO"].keys()
    for segment in meta["HALO"][flight]["segments"]
    if "circle" in segment["kinds"]
]
# %%
ds = ds.rename(
    {
        "latitude": "lat",
        "longitude": "lon",
    }
)

# %%
results_sonde = []
results_circle = []
for circle_name in [x for x in set(np.unique(ds.circle_id.values)) if x]:
    cds = ds.where(ds.circle_id == circle_name, drop=True)
    print(circle_name)
    circle_details = [
        segment for segment in details if segment["segment_id"] == circle_name
    ][0]
    flight_id = circle_name.split("_")[0]
    circle_results = []

    circle_results_sonde = []
    for platform in ["orc", "obs", "an", "fg"]:
        circle_ds = cds.sel(analysis_type=platform)
        circle_ds = circle_ds.assign(
            lat=circle_ds["lat"]
            .broadcast_like(circle_ds["p"])
            .transpose("profile", "p"),
            lon=circle_ds["lon"]
            .broadcast_like(circle_ds["p"])
            .transpose("profile", "p"),
        )

        for var in circle_ds.variables:
            circle_ds[var].attrs["standard_name"] = circle_ds[var].attrs.get(
                "standard_name", ""
            )
            circle_ds[var].attrs["long_name"] = circle_ds[var].attrs.get(
                "long_name", ""
            )

        circle = circles.Circle(
            circle_ds=circle_ds,
            clon=circle_details["clon"],
            clat=circle_details["clat"],
            crad=circle_details["radius"],
            flight_id=flight_id,
            platform_id=platform,
            segment_id=circle_name + platform,
            alt_dim="p",
            sonde_dim="profile",
        )

        circle = circle.get_xy_coords_for_circles()
        circle = circle.drop_vars()
        circle = circle.interpolate_na_sondes()
        circle = circle.extrapolate_na_sondes(max_alt=90000)
        circle = circle.apply_fit2d()
        circle = circle.add_divergence()
        circle = circle.add_vorticity()
        circle = circle.add_omega()

        res = circle.circle_ds.copy()
        sonde_vars = ["u", "v", "q", "ta", "p", "x", "y"]
        if platform != "orc":
            res["q"] = res["q"] / 1000
        circle_results_sonde.append(
            res[sonde_vars]
            .drop_vars("circle_id")
            .swap_dims({"profile": "launch_time"})
            .sortby("launch_time")
            .assign_coords(
                analysis_type=platform,
                circle_id_sonde=circle_name,
            )
        )
        circle_results.append(
            res.reset_coords()
            .drop_vars(
                sonde_vars + ["lat", "sonde_id", "launch_time", "lon", "circle_id"]
            )
            .assign_coords(
                circle_id=circle_name,
                analysis_type=platform,
            )
            .expand_dims("circle_id")
        )

    results_sonde.append(xr.concat(circle_results_sonde, dim="analysis_type"))
    results_circle.append(xr.concat(circle_results, dim="analysis_type"))


# %%
ds_sonde = xr.concat(circle_results_sonde, dim="analysis_type")
ds_circle = xr.concat(circle_results, dim="analysis_type")
# %%

# %%
circle_sondes = xr.combine_by_coords(
    results_sonde,
).swap_dims({"launch_time": "profile"})
circle_circle = xr.combine_by_coords(
    results_circle,
)
final = xr.merge(
    [circle_sondes, circle_circle],
)


# %%
def add_attrs(ds, var, attrs):
    ds[var].attrs.update(attrs)
    if var in ds.coords:
        print(var)
        ds = ds.assign_coords({var: (ds[var].dims, ds[var].values, attrs)})
    return ds


# %% add var attrubutes

for var, attrs in attr_dict.items():
    print(var)
    final = add_attrs(final, var, attrs)
for var, attrs in circle_attr_dict.items():
    final = add_attrs(final, var, attrs)
# %%

# %%

# %%
final.to_netcdf("ds_assimilated_omega.nc")
