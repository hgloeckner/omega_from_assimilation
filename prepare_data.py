# %%
import xarray as xr
import pydropsonde.circles as circles
import numpy as np
from orcestra import get_flight_segments
from pydropsonde.helper.xarray_helper import write_ds


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
ds = (
    ds.where(ds.launch_datetime != "--", drop=True)
    .rename({orig: new_var for orig, new_var in zip(original_vars, new_vars)})
    .dropna(dim="profile", how="all", subset=["u_obs", "v_obs"], thresh=400)
)

# %%
lev4 = xr.open_dataset(
    "ipfs://QmP73Kosem4exJcZXxG8vpN4YLqaepoZSPWwnQ9N1xffus/Level_4/PERCUSION_Level_4.zarr",
    engine="zarr",
)
lev3 = (
    xr.open_dataset(
        "ipfs://QmP73Kosem4exJcZXxG8vpN4YLqaepoZSPWwnQ9N1xffus/Level_3/PERCUSION_Level_3.zarr",
        engine="zarr",
    )
    .drop_vars(["bin_average_time"])
    .interpolate_na(dim="altitude")
)

# %%
ds = ds.assign(launch_time=ds.launch_time.astype(np.datetime64))

# %%

from scipy.stats import pearsonr


def get_correlations(ds, lev4, var="u_obs", l4_var="u"):
    corr = {}
    fail = []
    for profile in ds.profile.values:
        ds_profile = (
            ds.sel(profile=profile)
            .swap_dims({"level": "p_obs"})
            .dropna(dim="p_obs")
            .sortby("p_obs", ascending=False)
        )
        time = ds_profile.launch_time.values
        proxy_l4 = lev4.swap_dims({"sonde": "launch_time"}).sel(
            launch_time=slice(
                time - np.timedelta64(15, "m"), time + np.timedelta64(15, "m")
            )
        )
        corr[profile] = {}
        if proxy_l4.sizes["launch_time"] == 0:
            print(f"no proxy for {profile} at {time}")
            fail.append(profile)
            continue
        for launch_time in proxy_l4.launch_time:
            proxy_profile = (
                proxy_l4.sel(launch_time=launch_time)
                .swap_dims({"altitude": "p"})
                .dropna(dim="p", subset=["p"])
                .sortby("p", ascending=False)
            )
            if (proxy_profile.sizes["p"] > 0) and (ds_profile.sizes["p_obs"] > 0):
                proxy_profile = proxy_profile.sel(
                    p=ds_profile.p_obs.values, method="nearest"
                )
            else:
                print(f"no sonde for {profile} at {launch_time.values}")
                continue
            square_diff = np.mean(
                (
                    (ds_profile[var].values - proxy_profile[l4_var].values) ** 2
                    / (
                        (
                            np.abs(ds_profile[var].values)
                            + np.abs(proxy_profile[l4_var].values)
                        )
                        / 2
                    )
                )
            )

            corr[profile][str(proxy_profile.sonde_id.values)] = (
                pearsonr(ds_profile[var].values, proxy_profile[l4_var].values),
                square_diff,
            )
    return corr, fail


res = {}
fail = {}
for var in ["u", "v"]:
    res[var], fail[var] = get_correlations(ds, lev3, var=f"{var}_obs", l4_var=var)


# %%
def get_means(corr_dict):
    mean_corr = {}
    for profile in corr_dict["u"].keys():
        mean_corr[profile] = {}
        for sonde_id in corr_dict["u"][profile].keys():
            mean_corr[profile][sonde_id] = {
                "mean_corr": np.mean(
                    [
                        corr_dict["u"][profile][sonde_id][0][0],
                        corr_dict["v"][profile][sonde_id][0][0],
                    ]
                ),
                "mean_err": np.mean(
                    [
                        corr_dict["u"][profile][sonde_id][1],
                        corr_dict["v"][profile][sonde_id][1],
                    ]
                ),
            }

    return mean_corr


mean_corr = {}
mean_corr = get_means(res)


# %%
def get_best_match(mean_corr, match_var="mean_err"):
    best_match = {}
    for key, value in mean_corr.items():
        print(key, value)
        if value:
            best_sonde = min(value.items(), key=lambda x: x[1][match_var])
            best_match[key] = best_sonde[0]
        else:
            best_match[key] = ""
    return best_match


best_match = get_best_match(mean_corr, match_var="mean_err")


# %%
def find_ifs_sondes_for_duplicate(best_sondes_dict, sonde_id):
    profiles = []
    for profile, sonde in best_sondes_dict.items():
        if sonde == sonde_id:
            profiles.append(profile)
    return profiles


def find_best_matches_for_duplicates(
    mean_corr, duplicates, best_match, known_matches=None
):
    for duplicate_id in duplicates:
        if duplicate_id:
            profiles = find_ifs_sondes_for_duplicate(
                best_sondes_dict=best_match, sonde_id=duplicate_id
            )
            try:
                true_profile = profiles[
                    np.argmin(
                        [
                            mean_corr[profile][duplicate_id]["mean_err"]
                            for profile in profiles
                        ]
                    )
                ]
            except ValueError:
                print(profiles)
            for profile in profiles:
                if profile != true_profile:
                    profile_err = {}
                    for sonde_id in mean_corr[profile].keys():
                        if (sonde_id == duplicate_id) or (
                            np.isin(sonde_id, known_matches)
                        ):
                            profile_err[sonde_id] = mean_corr[profile][sonde_id][
                                "mean_err"
                            ]
                            mean_corr[profile][duplicate_id]["mean_err"] = np.inf
                    if np.all(np.isin(list(mean_corr[profile].keys()), known_matches)):
                        best_match[profile] = ""
                    else:
                        best_sonde = min(
                            mean_corr[profile].items(),
                            key=lambda x: x[1]["mean_err"],
                        )
                        best_match[profile] = best_sonde[0]
                    for sonde_id in profile_err.keys():
                        mean_corr[profile][sonde_id]["mean_err"] = profile_err[sonde_id]
    return best_match


# %%
for i in range(20):
    duplicates = {
        x for x in list(best_match.values()) if list(best_match.values()).count(x) > 1
    }
    known_matches = np.unique(list(best_match.values()))
    new_best = find_best_matches_for_duplicates(
        mean_corr, duplicates, best_match, known_matches=known_matches
    )


# %%
assert ds.sizes["profile"] == len(new_best)
assert sorted(list(best_match.keys())) == list(best_match.keys())

ds = ds.assign(sonde_id=("profile", list(best_match.values())))
ds = ds.where(ds.sonde_id, drop=True)

# %%


def get_circle_id_from_sonde_id(lev4, sonde_id):
    try:
        sonde_idx = (
            lev4.swap_dims({"sonde": "sonde_id"})
            .get_index("sonde_id")
            .get_loc(sonde_id)
        )
    except KeyError:
        return ""
    circle_end = lev4.sondes_per_circle.cumsum()
    cidx = np.argmin(np.abs(sonde_idx - circle_end.values))
    if circle_end[cidx] > sonde_idx:
        return cidx
    else:
        return cidx + 1


# %% add orcestra data

orc_ifs = (
    lev4[["u", "v", "q", "ta", "p", "sonde_id"]]
    .rename(
        {
            "u": "u_orc",
            "v": "v_orc",
            "q": "q_orc",
            "ta": "ta_orc",
            "p": "p_orc",
        }
    )
    .swap_dims({"sonde": "sonde_id"})
)
ds = ds.swap_dims({"profile": "sonde_id"})
# %%
sonde_ids_in_both = set(orc_ifs.sonde_id.values).intersection(set(ds.sonde_id.values))
orc_ifs = orc_ifs.sel(sonde_id=list(sonde_ids_in_both))
# %%
full_ds = xr.merge(
    [
        ds.drop_vars(["launch_time"]),
        orc_ifs.swap_dims({"altitude": "level"})
        .interp(level=ds.level)
        .drop_vars(["level"]),
    ]
)
# %%  make p dimension
from xhistogram.xarray import histogram


def make_p_dimension(ds, analysis_type):
    vars = [var for var in ds.variables if var.endswith(analysis_type)]
    p_var = f"p_{analysis_type}"
    p_bins = np.logspace(np.log(102000), np.log(14000), num=10000, base=np.e)[::-1]

    mean_ds = []
    count_dict = {}
    for var in vars:
        if var in ds.variables and var not in ds.dims:
            count_dict[var] = histogram(
                ds[p_var].where(~np.isnan(ds[var])),
                bins=p_bins,
                dim=["level"],
            )

            new_ds = (
                histogram(
                    ds[p_var].where(~np.isnan(ds[var])),
                    bins=p_bins,
                    dim=["level"],
                    weights=ds[var].astype(np.float64).where(~np.isnan(ds[var])),
                )  # casting necessary for time
                / count_dict[var]
            )
            new_ds.name = var
            mean_ds.append(new_ds.copy())

    binned = xr.merge(mean_ds).rename({f"{p_var}_bin": "p"})
    binned = binned.assign(p=("p", np.log(binned.p.values)))
    binned = binned.interpolate_na(dim="p")
    binned = binned.assign(p=("p", np.exp(binned.p.values)))
    for var in binned.variables:
        if var in ds.variables:
            binned[var].attrs = ds[var].attrs
    return binned


res = {}
for analysis_type in ["obs", "an", "fg", "orc"]:
    res[analysis_type] = []
    res[analysis_type].append(
        make_p_dimension(full_ds, analysis_type)
        .rename({f"{var}_{analysis_type}": var for var in ["u", "v", "q", "ta"]})
        .assign(analysis_type=analysis_type)
        .drop_vars(f"p_{analysis_type}")
    )
# %%
result_ds = xr.concat(
    [xr.merge(res[analysis_type]) for analysis_type in ["obs", "an", "fg", "orc"]],
    dim="analysis_type",
).dropna(dim="p", how="all")


# %%
def fix_q_unit(ds):
    for var in ds.variables:
        if ds[var].attrs.get("units") == "g/kg":
            print(var)
            ds[var] = ds[var] / 1000  # convert g kg-1 to kg kg-1
            ds[var].attrs["units"] = "kg kg-1"
    return ds


circle_ids = [
    lev4.circle_id.sel(circle=get_circle_id_from_sonde_id(lev4, sonde_id)).values
    for sonde_id in result_ds.sonde_id.values
]
result_ds = result_ds.assign_coords(circle_id=("sonde_id", circle_ids))
# %%
#
write_ds(
    result_ds.swap_dims({"sonde_id": "profile"}),
    dir=".",
    filename="ds_assimilated_and_orcestra.nc",
    object_dims=("profile", "analysis_type"),
    alt_dim="p",
)
# %%
ds = xr.open_dataset("ds_assimilated_and_orcestra.nc")

# %%
profile = 20
result_ds.isel(sonde_id=profile).sel(analysis_type="orc").v.plot()

result_ds.isel(sonde_id=profile).sel(analysis_type="obs").v.plot()

# %%
