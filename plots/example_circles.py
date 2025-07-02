# %%
import matplotlib.pyplot as plt
import xarray as xr
import seaborn as sns


ds = xr.open_dataset("../ds_assimilated_omega.nc")
lev4 = xr.open_dataset(
    "ipfs://QmP73Kosem4exJcZXxG8vpN4YLqaepoZSPWwnQ9N1xffus/Level_4/PERCUSION_Level_4.zarr",
    engine="zarr",
).swap_dims({"circle": "circle_id"})
# %%
circles = ["HALO-20240831a_bb86", "HALO-20240924a_819e"]
var = "omega"
fig, axes = plt.subplots(ncols=2, figsize=(6, 3))

for ax, circle in zip(axes, circles):
    for analysis_type in ds.analysis_type.values:
        ds_circle = ds.sel(circle_id=circle, analysis_type=analysis_type).assign(
            p=ds.sel(circle_id=circle, analysis_type=analysis_type).p / 100
        )
        if circle == circles[1]:
            ds_circle[var].plot(
                y="p",
                ax=ax,
                label=analysis_type,
            )
        else:
            ds_circle[var].plot(
                y="p",
                ax=ax,
            )
    ax.scatter(
        lev4.sel(circle_id=circle)[var].values,
        lev4.sel(circle_id=circle).p_mean.values / 100,
        color="k",
        marker="x",
        label="PERCUSION",
        s=2,
    )
    ax.set_title(circle)
    ax.invert_yaxis()
    ax.set_xlabel("Omega / hPa hr-1")
    ax.set_ylabel("")
axes[0].set_ylabel("Pressure / hPa")
axes[1].legend()
axes[1].set_yticklabels([])
sns.despine(offset=10)
fig.savefig("ifs_compare.pdf", bbox_inches="tight")
